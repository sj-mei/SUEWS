from typing import Optional
from pydantic import ConfigDict, BaseModel, Field, model_validator
from .type import RefValue, Reference, FlexibleRefValue
from .profile import HourlyProfile
from .type import init_df_state
from .surface import (
    SurfaceType,
    SurfaceProperties,
    PavedProperties,
    BldgsProperties,
    BsoilProperties,
    WaterProperties,
    VerticalLayers,
)
from .human_activity import AnthropogenicEmissions, IrrigationParams
from .hydro import (
    WaterDistribution,
)
from .state import InitialStates

import pandas as pd
from typing import List, Literal, Union, Dict, Tuple

from datetime import datetime
from timezonefinder import TimezoneFinder
from pytz import timezone
from pytz.exceptions import AmbiguousTimeError, NonExistentTimeError
import pytz
import warnings


class VegetationParams(BaseModel):
    porosity_id: FlexibleRefValue(int) = Field(
        description="Initial porosity for deciduous trees", 
        json_schema_extra={"unit": "dimensionless"}
    )
    gdd_id: FlexibleRefValue(int) = Field(
        description="Growing degree days ID", 
        json_schema_extra={"unit": "degC d"}
    )
    sdd_id: FlexibleRefValue(int) = Field(
        description="Senescence degree days ID", 
        json_schema_extra={"unit": "degC d"}
    )
    lai: Dict[str, Union[FlexibleRefValue(float), List[FlexibleRefValue(float)]]] = Field(
        description="Leaf area index parameters",
        json_schema_extra={"unit": "m^2 m^-2"}
    )
    ie_a: FlexibleRefValue(float) = Field(
        description="Irrigation efficiency coefficient a", 
        json_schema_extra={"unit": "dimensionless"}
    )
    ie_m: FlexibleRefValue(float) = Field(
        description="Irrigation efficiency coefficient m", 
        json_schema_extra={"unit": "dimensionless"}
    )

    ref: Optional[Reference] = None


class Conductance(BaseModel):
    g_max: FlexibleRefValue(float) = Field(
        default=40.0,
        description="Maximum surface conductance for photosynthesis", json_schema_extra={"unit": "mm s^-1"}
    )
    g_k: FlexibleRefValue(float) = Field(
        default=0.6,
        description="Conductance parameter related to incoming solar radiation", json_schema_extra={"unit": "dimensionless"}
    )
    g_q_base: FlexibleRefValue(float) = Field(
        default=0.03,
        description="Base value for conductance parameter related to vapor pressure deficit", json_schema_extra={"unit": "kPa^-1"}
    )
    g_q_shape: FlexibleRefValue(float) = Field(
        default=0.9,
        description="Shape parameter for conductance related to vapor pressure deficit", json_schema_extra={"unit": "dimensionless"}
    )
    g_t: FlexibleRefValue(float) = Field(
        default=30.0,
        description="Conductance parameter related to air temperature", json_schema_extra={"unit": "degC"}
    )
    g_sm: FlexibleRefValue(float) = Field(
        default=0.5,
        description="Conductance parameter related to soil moisture", json_schema_extra={"unit": "dimensionless"}
    )
    kmax: FlexibleRefValue(float) = Field(
        default=1200.0,
        description="Maximum incoming shortwave radiation", json_schema_extra={"unit": "W m^-2"}
    )
    s1: FlexibleRefValue(float) = Field(
        default=0.2,
        description="Lower soil moisture threshold for conductance response", 
        json_schema_extra={"unit": "dimensionless"}
    )
    s2: FlexibleRefValue(float) = Field(
        default=0.5,
        description="Parameter related to soil moisture dependence", json_schema_extra={"unit": "mm"}
    )
    tl: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Lower air temperature threshold for conductance response", json_schema_extra={"unit": "degC"}
    )
    th: FlexibleRefValue(float) = Field(
        default=50.0,
        description="Upper air temperature threshold for conductance response", json_schema_extra={"unit": "degC"}
    )

    ref: Optional[Reference] = Reference(ref="Test ref", DOI="test doi", ID="test id")

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """
        Convert conductance parameters to DataFrame state format.

        Args:
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            pd.DataFrame: DataFrame containing conductance parameters.
        """

        df_state = init_df_state(grid_id)

        scalar_params = {
            "g_max": self.g_max,
            "g_k": self.g_k,
            "g_q_base": self.g_q_base,
            "g_q_shape": self.g_q_shape,
            "g_t": self.g_t,
            "g_sm": self.g_sm,
            "kmax": self.kmax,
            "s1": self.s1,
            "s2": self.s2,
            "tl": self.tl,
            "th": self.th,
        }

        for param_name, value in scalar_params.items():
            val = value.value if isinstance(value, RefValue) else value
            df_state.loc[grid_id, (param_name, "0")] = val

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "Conductance":
        """
        Reconstruct Conductance from a DataFrame state format.

        Args:
            df: DataFrame containing conductance parameters
            grid_id: Grid ID for the DataFrame index

        Returns:
            Conductance: Instance of Conductance
        """
        scalar_params = {
            "g_max": df.loc[grid_id, ("g_max", "0")],
            "g_k": df.loc[grid_id, ("g_k", "0")],
            "g_q_base": df.loc[grid_id, ("g_q_base", "0")],
            "g_q_shape": df.loc[grid_id, ("g_q_shape", "0")],
            "g_t": df.loc[grid_id, ("g_t", "0")],
            "g_sm": df.loc[grid_id, ("g_sm", "0")],
            "kmax": df.loc[grid_id, ("kmax", "0")],
            "s1": df.loc[grid_id, ("s1", "0")],
            "s2": df.loc[grid_id, ("s2", "0")],
            "tl": df.loc[grid_id, ("tl", "0")],
            "th": df.loc[grid_id, ("th", "0")],
        }

        # Convert scalar parameters to RefValue
        scalar_params = {
            key: RefValue(value) for key, value in scalar_params.items()
        }

        return cls(**scalar_params)


class LAIPowerCoefficients(BaseModel):
    growth_lai: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Power coefficient for LAI in growth equation (LAIPower[1])",
        json_schema_extra={"unit": "dimensionless"}
    )
    growth_gdd: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Power coefficient for GDD in growth equation (LAIPower[2])",
        json_schema_extra={"unit": "dimensionless"}
    )
    senescence_lai: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Power coefficient for LAI in senescence equation (LAIPower[3])",
        json_schema_extra={"unit": "dimensionless"}
    )
    senescence_sdd: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Power coefficient for SDD in senescence equation (LAIPower[4])",
        json_schema_extra={"unit": "dimensionless"}
    )

    ref: Optional[Reference] = None

    def to_list(self) -> List[float]:
        """Convert to list format for Fortran interface"""
        return [
            self.growth_lai,
            self.growth_gdd,
            self.senescence_lai,
            self.senescence_sdd,
        ]

    def to_df_state(self, grid_id: int, veg_idx: int) -> pd.DataFrame:
        """Convert LAI power coefficients to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            veg_idx: Vegetation index (0: EVETR, 1: DECTR, 2: GRASS)

        Returns:
            pd.DataFrame: DataFrame containing LAI power coefficients
        """
        df_state = init_df_state(grid_id)

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, indices: Tuple, value: float):
            idx_str = str(indices)
            if (col_name, idx_str) not in df_state.columns:
                # df_state[(col_name, idx_str)] = np.nan
                df_state[(col_name, idx_str)] = None
            df_state.at[grid_id, (col_name, idx_str)] = value

        # Set power coefficients in order
        for i, value in enumerate(self.to_list()):
            val = value.value if isinstance(value, RefValue) else value
            set_df_value("laipower", (i, veg_idx), val)

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, veg_idx: int
    ) -> "LAIPowerCoefficients":
        """
        Reconstruct LAIPowerCoefficients from DataFrame state format.

        Args:
            df: DataFrame containing LAI power coefficients
            grid_id: Grid ID for the DataFrame index
            veg_idx: Vegetation index (0: EVETR, 1: DECTR, 2: GRASS)

        Returns:
            LAIPowerCoefficients: Instance of LAIPowerCoefficients
        """
        # Map each coefficient to its corresponding index
        coefficients = [
            RefValue(df.loc[grid_id, ("laipower", f"(0, {veg_idx})")]),
            RefValue(df.loc[grid_id, ("laipower", f"(1, {veg_idx})")]),
            RefValue(df.loc[grid_id, ("laipower", f"(2, {veg_idx})")]),
            RefValue(df.loc[grid_id, ("laipower", f"(3, {veg_idx})")]),
        ]

        # Return the instance with coefficients
        return cls(
            growth_lai=coefficients[0],
            growth_gdd=coefficients[1],
            senescence_lai=coefficients[2],
            senescence_sdd=coefficients[3],
        )


class LAIParams(BaseModel):
    baset: FlexibleRefValue(float) = Field(
        default=10.0,
        description="Base temperature for initiating growing degree days (GDD) for leaf growth",
        json_schema_extra={"unit": "degC"}
    )
    gddfull: FlexibleRefValue(float) = Field(
        default=100.0,
        description="Growing degree days (GDD) needed for full capacity of LAI",
        json_schema_extra={"unit": "degC*day"}
    )
    basete: FlexibleRefValue(float) = Field(
        default=10.0,
        description="Base temperature for initiating senescence degree days (SDD) for leaf off",
        json_schema_extra={"unit": "degC"}
    )
    sddfull: FlexibleRefValue(float) = Field(
        default=100.0,
        description="Senescence degree days (SDD) needed to initiate leaf off",
        json_schema_extra={"unit": "degC*day"}
    )
    laimin: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Leaf-off wintertime LAI value", json_schema_extra={"unit": "m^2 m^-2"}
    )
    laimax: FlexibleRefValue(float) = Field(
        default=10.0,
        description="Full leaf-on summertime LAI value", json_schema_extra={"unit": "m^2 m^-2"}
    )
    laipower: LAIPowerCoefficients = Field(
        default_factory=LAIPowerCoefficients,
        description="LAI calculation power parameters for growth and senescence",
    )
    laitype: FlexibleRefValue(int) = Field(
        default=0,
        description="LAI calculation choice (0: original, 1: new high latitude)",
        json_schema_extra={"unit": "dimensionless"}
    )

    ref: Optional[Reference] = None

    @model_validator(mode="after")
    def validate_lai_ranges(self) -> "LAIParams":
        laimin_val = self.laimin.value if isinstance(self.laimin, RefValue) else self.laimin
        laimax_val = self.laimax.value if isinstance(self.laimax, RefValue) else self.laimax
        baset_val = self.baset.value if isinstance(self.baset, RefValue) else self.baset
        gddfull_val = self.gddfull.value if isinstance(self.gddfull, RefValue) else self.gddfull
        
        if laimin_val > laimax_val:
            raise ValueError(
                f"laimin ({laimin_val})must be less than or equal to laimax ({laimax_val})."
            )
        if baset_val > gddfull_val:
            raise ValueError(
                f"baset {baset_val} must be less than gddfull ({gddfull_val})."
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

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, indices: Union[Tuple, int], value: float):
            idx_str = str(indices) if isinstance(indices, int) else str(indices)
            if (col_name, idx_str) not in df_state.columns:
                # df_state[(col_name, idx_str)] = np.nan
                df_state[(col_name, idx_str)] = None
            df_state.at[grid_id, (col_name, idx_str)] = value

        # Set basic LAI parameters
        lai_params = {
            "baset": self.baset,
            "gddfull": self.gddfull,
            "basete": self.basete,
            "sddfull": self.sddfull,
            "laimin": self.laimin,
            "laimax": self.laimax,
            "laitype": self.laitype,
        }

        for param, value in lai_params.items():
            val = value.value if isinstance(value, RefValue) else value
            set_df_value(param, (veg_idx,), val)

        # Add LAI power coefficients using the LAIPowerCoefficients to_df_state method
        if self.laipower:
            df_power = self.laipower.to_df_state(grid_id, veg_idx)
            # Merge power coefficients into main DataFrame
            for col in df_power.columns:
                if col[0] != "grid_iv":  # Skip the grid_iv column
                    df_state[col] = df_power[col]

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int
    ) -> "LAIParams":
        """
        Reconstruct LAIParams from DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing LAI parameters.
            grid_id (int): Grid ID for the DataFrame index.
            surf_idx (int): Surface index for vegetation (2: EVETR, 3: DECTR, 4: GRASS).

        Returns:
            LAIParams: Instance of LAIParams.
        """
        # Adjust index for vegetation surfaces (surface index - 2)
        veg_idx = surf_idx - 2

        # Helper function to extract values from DataFrame
        def get_df_value(col_name: str, indices: Union[Tuple, int]) -> float:
            idx_str = str(indices) if isinstance(indices, int) else str(indices)
            return df.loc[grid_id, (col_name, idx_str)]

        # Extract basic LAI parameters
        lai_params = {
            "baset": get_df_value("baset", (veg_idx,)),
            "gddfull": get_df_value("gddfull", (veg_idx,)),
            "basete": get_df_value("basete", (veg_idx,)),
            "sddfull": get_df_value("sddfull", (veg_idx,)),
            "laimin": get_df_value("laimin", (veg_idx,)),
            "laimax": get_df_value("laimax", (veg_idx,)),
            "laitype": int(get_df_value("laitype", (veg_idx,))),
        }

        # Convert scalar parameters to RefValue
        lai_params = {key: RefValue(value) for key, value in lai_params.items()}

        # Extract LAI power coefficients
        laipower = LAIPowerCoefficients.from_df_state(df, grid_id, veg_idx)

        return cls(**lai_params, laipower=laipower)


class VegetatedSurfaceProperties(SurfaceProperties):
    alb: FlexibleRefValue(float) = Field(
        ge=0, le=1,
        description="Albedo", json_schema_extra={"unit": "dimensionless"},
        default=0.2
    )
    alb_min: FlexibleRefValue(float) = Field(
        ge=0, le=1,
        description="Minimum albedo", json_schema_extra={"unit": "dimensionless"},
        default=0.2
    )
    alb_max: FlexibleRefValue(float) = Field(
        ge=0, le=1,
        description="Maximum albedo", json_schema_extra={"unit": "dimensionless"},
        default=0.3
    )
    beta_bioco2: FlexibleRefValue(float) = Field(
        default=0.6, description="Biogenic CO2 exchange coefficient", json_schema_extra={"unit": "dimensionless"}
    )
    beta_enh_bioco2: FlexibleRefValue(float) = Field(
        default=0.7,
        description="Enhanced biogenic CO2 exchange coefficient", json_schema_extra={"unit": "dimensionless"},
    )
    alpha_bioco2: FlexibleRefValue(float) = Field(
        default=0.8, description="Biogenic CO2 exchange coefficient", json_schema_extra={"unit": "dimensionless"}
    )
    alpha_enh_bioco2: FlexibleRefValue(float) = Field(
        default=0.9,
        description="Enhanced biogenic CO2 exchange coefficient", json_schema_extra={"unit": "dimensionless"},
    )
    resp_a: FlexibleRefValue(float) = Field(
        default=1.0, description="Respiration coefficient", json_schema_extra={"unit": "umol m^-2 s^-1"}
    )
    resp_b: FlexibleRefValue(float) = Field(
        default=1.1, description="Respiration coefficient", json_schema_extra={"unit": "dimensionless"}
    )
    theta_bioco2: FlexibleRefValue(float) = Field(
        default=1.2, description="Biogenic CO2 exchange coefficient", json_schema_extra={"unit": "dimensionless"}
    )
    maxconductance: FlexibleRefValue(float) = Field(
        default=0.5, description="Maximum surface conductance", json_schema_extra={"unit": "mm s^-1"}
    )
    min_res_bioco2: FlexibleRefValue(float) = Field(
        default=0.1, description="Minimum respiratory biogenic CO2", json_schema_extra={"unit": "umol m^-2 s^-1"}
    )
    lai: LAIParams = Field(
        default_factory=LAIParams, description="Leaf area index parameters"
    )
    ie_a: FlexibleRefValue(float) = Field(
        default=0.5,
        description="Irrigation efficiency coefficient-automatic", json_schema_extra={"unit": "dimensionless"},
    )
    ie_m: FlexibleRefValue(float) = Field(
        default=0.6,
        description="Irrigation efficiency coefficient-manual", json_schema_extra={"unit": "dimensionless"},
    )

    ref: Optional[Reference] = None

    @model_validator(mode="after")
    def validate_albedo_range(self) -> "VegetatedSurfaceProperties":
        alb_min_val = self.alb_min.value if isinstance(self.alb_min, RefValue) else self.alb_min
        alb_max_val = self.alb_max.value if isinstance(self.alb_max, RefValue) else self.alb_max
        
        if alb_min_val > alb_max_val:
            raise ValueError(
                f"alb_min (input {alb_min_val}) must be less than or equal to alb_max (entered {alb_max_val})."
            )
        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert vegetated surface properties to DataFrame state format."""
        # Get base properties
        df_state = super().to_df_state(grid_id)

        # Add vegetation-specific properties
        surf_idx = self.get_surface_index()

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, idx_str: str, value: float):
            if (col_name, idx_str) not in df_state.columns:
                # df_state[(col_name, idx_str)] = np.nan
                df_state[(col_name, idx_str)] = None
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # add ordinary float properties
        for attr in [
            "alb",
            # "alb_min",
            # "alb_max",
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
            field_val = getattr(self, attr)
            val = field_val.value if isinstance(field_val, RefValue) else field_val
            set_df_value(attr, f"({surf_idx-2},)", val)

        df_lai = self.lai.to_df_state(grid_id, surf_idx)
        df_state = pd.concat([df_state, df_lai], axis=1).sort_index(axis=1)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int, surf_idx: int) -> "VegetatedSurfaceProperties":
        """Reconstruct vegetated surface properties from DataFrame state format."""
        instance = super().from_df_state(df, grid_id, surf_idx)
        # add ordinary float properties
        for attr in [
            "alb",
            # "alb_min",
            # "alb_max",
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
            setattr(instance, attr, RefValue(df.loc[grid_id, (attr, f"({surf_idx-2},)")]))

        instance.lai = LAIParams.from_df_state(df, grid_id, surf_idx)

        return instance


class EvetrProperties(VegetatedSurfaceProperties):  # TODO: Move waterdist VWD here?
    alb: FlexibleRefValue(float) = Field(
        ge=0, le=1,
        default=0.2,
        description="Albedo", json_schema_extra={"unit": "dimensionless"}
    )
    faievetree: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Frontal area index of evergreen trees", json_schema_extra={"unit": "dimensionless"}
    )
    evetreeh: FlexibleRefValue(float) = Field(
        default=15.0,
        description="Evergreen tree height", json_schema_extra={"unit": "m"}
    )
    _surface_type: Literal[SurfaceType.EVETR] = SurfaceType.EVETR
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.EVETR),
        description="Water distribution for evergreen trees",
    )

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert evergreen tree properties to DataFrame state format."""
        # Get base properties from parent
        df_state = super().to_df_state(grid_id)
        surf_idx = self.get_surface_index()

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = f"({surf_idx},)"
            if (col_name, idx_str) not in df_state.columns:
                # df_state[(col_name, idx_str)] = np.nan
                df_state[(col_name, idx_str)] = None
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # Add all non-inherited properties
        list_properties = ["faievetree", "evetreeh"]
        for attr in list_properties:
            field_val = getattr(self, attr)
            val = field_val.value if isinstance(field_val, RefValue) else field_val
            df_state.loc[grid_id, (attr, "0")] = val

        # specific properties
        df_state.loc[grid_id, ("alb", "(2,)")] = self.alb.value if isinstance(self.alb, RefValue) else self.alb
        df_state.loc[grid_id, ("albmin_evetr", "0")] = self.alb_min.value if isinstance(self.alb_min, RefValue) else self.alb_min
        df_state.loc[grid_id, ("albmax_evetr", "0")] = self.alb_max.value if isinstance(self.alb_max, RefValue) else self.alb_max

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "EvetrProperties":
        """Reconstruct evergreen tree properties from DataFrame state format."""
        surf_idx = 2
        instance = super().from_df_state(df, grid_id, surf_idx)

        instance.alb = RefValue(df.loc[grid_id, ("alb", "(2,)")])
        instance.faievetree = RefValue(df.loc[grid_id, ("faievetree", "0")])
        instance.evetreeh = RefValue(df.loc[grid_id, ("evetreeh", "0")])

        instance.alb_min = RefValue(df.loc[grid_id, ("albmin_evetr", "0")])
        instance.alb_max = RefValue(df.loc[grid_id, ("albmax_evetr", "0")])

        return instance


class DectrProperties(VegetatedSurfaceProperties):
    alb: FlexibleRefValue(float) = Field(
        ge=0, le=1,
        default=0.2,
        description="Albedo", json_schema_extra={"unit": "dimensionless"}
    )
    faidectree: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Frontal area index of deciduous trees", json_schema_extra={"unit": "dimensionless"}
    )
    dectreeh: FlexibleRefValue(float) = Field(
        default=15.0,
        description="Deciduous tree height", json_schema_extra={"unit": "m"}
    )
    pormin_dec: FlexibleRefValue(float) = Field(
        ge=0.1, le=0.9,
        default=0.2,
        description="Minimum porosity", json_schema_extra={"unit": "dimensionless"}
    )  # pormin_dec cannot be less than 0.1 and greater than 0.9
    pormax_dec: FlexibleRefValue(float) = Field(
        ge=0.1, le=0.9,
        default=0.6,
        description="Maximum porosity", json_schema_extra={"unit": "dimensionless"}
    )  # pormax_dec cannot be less than 0.1 and greater than 0.9
    capmax_dec: FlexibleRefValue(float) = Field(
        default=100.0,
        description="Maximum water capacity", json_schema_extra={"unit": "mm"}
    )
    capmin_dec: FlexibleRefValue(float) = Field(
        default=10.0,
        description="Minimum water capacity", json_schema_extra={"unit": "mm"}
    )
    _surface_type: Literal[SurfaceType.DECTR] = SurfaceType.DECTR
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.DECTR),
        description="Water distribution for deciduous trees",
    )

    ref: Optional[Reference] = None

    @model_validator(mode="after")
    def validate_porosity_range(self) -> "DectrProperties":
        pormin_dec_val = self.pormin_dec.value if isinstance(self.pormin_dec, RefValue) else self.pormin_dec
        pormax_dec_val = self.pormax_dec.value if isinstance(self.pormax_dec, RefValue) else self.pormax_dec
        
        if pormin_dec_val >= pormax_dec_val:
            raise ValueError(
                f"pormin_dec ({pormin_dec_val}) must be less than pormax_dec ({pormax_dec_val})."
            )
        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert deciduous tree properties to DataFrame state format."""
        # Get base properties from parent
        df_state = super().to_df_state(grid_id)

        list_properties = [
            "faidectree",
            "dectreeh",
            "pormin_dec",
            "pormax_dec",
            "capmax_dec",
            "capmin_dec",
        ]
        # Add all non-inherited properties
        for attr in list_properties:
            field_val = getattr(self, attr)
            val = field_val.value if isinstance(field_val, RefValue) else field_val
            df_state.loc[grid_id, (attr, "0")] = val

        # specific properties
        df_state.loc[grid_id, ("alb", "(3,)")] = self.alb.value if isinstance(self.alb, RefValue) else self.alb
        df_state.loc[grid_id, ("albmin_dectr", "0")] = self.alb_min.value if isinstance(self.alb_min, RefValue) else self.alb_min
        df_state.loc[grid_id, ("albmax_dectr", "0")] = self.alb_max.value if isinstance(self.alb_max, RefValue) else self.alb_max

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "DectrProperties":
        """Reconstruct deciduous tree properties from DataFrame state format."""
        surf_idx = 3
        instance = super().from_df_state(df, grid_id, surf_idx)

        instance.alb = RefValue(df.loc[grid_id, ("alb", "(3,)")])
        instance.faidectree = RefValue(df.loc[grid_id, ("faidectree", "0")])
        instance.dectreeh = RefValue(df.loc[grid_id, ("dectreeh", "0")])
        instance.pormin_dec = RefValue(df.loc[grid_id, ("pormin_dec", "0")])
        instance.pormax_dec = RefValue(df.loc[grid_id, ("pormax_dec", "0")])
        instance.capmax_dec = RefValue(df.loc[grid_id, ("capmax_dec", "0")])
        instance.capmin_dec = RefValue(df.loc[grid_id, ("capmin_dec", "0")])

        instance.alb_min = RefValue(df.loc[grid_id, ("albmin_dectr", "0")])
        instance.alb_max = RefValue(df.loc[grid_id, ("albmax_dectr", "0")])

        return instance


class GrassProperties(VegetatedSurfaceProperties):
    alb: FlexibleRefValue(float) = Field(
        ge=0, le=1, default=0.2, description="Minimum albedo", json_schema_extra={"unit": "dimensionless"}
    )
    _surface_type: Literal[SurfaceType.GRASS] = SurfaceType.GRASS
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.GRASS),
        description="Water distribution for grass",
    )

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert grass properties to DataFrame state format."""
        # Get base properties from parent
        df_state = super().to_df_state(grid_id)

        # add specific properties
        df_state.loc[grid_id, ("alb", "(4,)")] = self.alb.value if isinstance(self.alb, RefValue) else self.alb
        df_state[("albmin_grass", "0")] = self.alb_min.value if isinstance(self.alb_min, RefValue) else self.alb_min
        df_state[("albmax_grass", "0")] = self.alb_max.value if isinstance(self.alb_max, RefValue) else self.alb_max

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "GrassProperties":
        """Reconstruct grass properties from DataFrame state format."""
        surf_idx = 4
        instance = super().from_df_state(df, grid_id, surf_idx)

        instance.alb = RefValue(df.loc[grid_id, ("alb", "(4,)")])
        instance.alb_min = RefValue(df.loc[grid_id, ("albmin_grass", "0")])
        instance.alb_max = RefValue(df.loc[grid_id, ("albmax_grass", "0")])

        return instance


class SnowParams(BaseModel):
    crwmax: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Maximum water holding capacity of snow", json_schema_extra={"unit": "mm"}
    )
    crwmin: FlexibleRefValue(float) = Field(
        default=0.05,
        description="Minimum water holding capacity of snow", json_schema_extra={"unit": "mm"}
    )
    narp_emis_snow: FlexibleRefValue(float) = Field(
        default=0.99,
        description="Snow surface emissivity", json_schema_extra={"unit": "dimensionless"}
    )
    preciplimit: FlexibleRefValue(float) = Field(
        default=2.2,
        description="Temperature threshold for snow vs rain precipitation", json_schema_extra={"unit": "degC"}
    )
    preciplimitalb: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Precipitation threshold for snow albedo aging", json_schema_extra={"unit": "mm"}
    )
    snowalbmax: FlexibleRefValue(float) = Field(
        default=0.85,
        description="Maximum snow albedo", json_schema_extra={"unit": "dimensionless"}
    )
    snowalbmin: FlexibleRefValue(float) = Field(
        default=0.4,
        description="Minimum snow albedo", json_schema_extra={"unit": "dimensionless"}
    )
    snowdensmin: FlexibleRefValue(float) = Field(
        default=100.0,
        description="Minimum snow density", json_schema_extra={"unit": "kg m^-3"}
    )
    snowdensmax: FlexibleRefValue(float) = Field(
        default=400.0,
        description="Maximum snow density", json_schema_extra={"unit": "kg m^-3"}
    )
    snowlimbldg: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Maximum snow depth limit on buildings", json_schema_extra={"unit": "m"}
    )
    snowlimpaved: FlexibleRefValue(float) = Field(
        default=0.1,
        description="Maximum snow depth limit on paved surfaces", json_schema_extra={"unit": "m"}
    )
    snowprof_24hr: HourlyProfile = Field(
        default_factory=HourlyProfile, description="24-hour snow profile"
    )
    tau_a: FlexibleRefValue(float) = Field(
        default=0.018,
        description="Time constant for snow albedo aging in cold snow", json_schema_extra={"unit": "dimensionless"}
    )
    tau_f: FlexibleRefValue(float) = Field(
        default=0.11,
        description="Time constant for snow albedo aging in melting snow", json_schema_extra={"unit": "dimensionless"}
    )
    tau_r: FlexibleRefValue(float) = Field(
        default=0.05,
        description="Time constant for snow albedo aging in refreezing snow", json_schema_extra={"unit": "dimensionless"}
    )
    tempmeltfact: FlexibleRefValue(float) = Field(
        default=0.12,
        description="Hourly temperature melt factor of snow", json_schema_extra={"unit": "mm K^-1 h^-1"}
    )
    radmeltfact: FlexibleRefValue(float) = Field(
        default=0.0016,
        description="Hourly radiation melt factor of snow", json_schema_extra={"unit": "mm W^-1 m^2 h^-1"}
    )

    ref: Optional[Reference] = None

    @model_validator(mode="after")
    def validate_all(self) -> "SnowParams":
        """
        Aggregate all validation checks for SnowParams,
        so multiple errors (if any) can be raised together
        """
        errors = []
        crwmin_val = self.crwmin.value if isinstance(self.crwmin, RefValue) else self.crwmin
        crwmax_val = self.crwmax.value if isinstance(self.crwmax, RefValue) else self.crwmax
        snowalbmin_val = self.snowalbmin.value if isinstance(self.snowalbmin, RefValue) else self.snowalbmin
        snowalbmax_val = self.snowalbmax.value if isinstance(self.snowalbmax, RefValue) else self.snowalbmax
        
        if crwmin_val >= crwmax_val:
            errors.append(
                f"crwmin ({crwmin_val}) must be less than crwmax ({crwmax_val})."
            )
        if snowalbmin_val >= snowalbmax_val:
            errors.append(
                f"snowalbmin ({snowalbmin_val}) must be less than snowalbmax ({snowalbmax_val})."
            )
        if errors:
            raise ValueError("\n".join(errors))

        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """
        Convert snow parameters to DataFrame state format.

        Args:
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            pd.DataFrame: DataFrame containing snow parameters.
        """

        df_state = init_df_state(grid_id)

        scalar_params = {
            "crwmax": self.crwmax,
            "crwmin": self.crwmin,
            "narp_emis_snow": self.narp_emis_snow,
            "preciplimit": self.preciplimit,
            "preciplimitalb": self.preciplimitalb,
            "snowalbmax": self.snowalbmax,
            "snowalbmin": self.snowalbmin,
            "snowdensmin": self.snowdensmin,
            "snowdensmax": self.snowdensmax,
            "snowlimbldg": self.snowlimbldg,
            "snowlimpaved": self.snowlimpaved,
            "tau_a": self.tau_a,
            "tau_f": self.tau_f,
            "tau_r": self.tau_r,
            "tempmeltfact": self.tempmeltfact,
            "radmeltfact": self.radmeltfact,
        }
        for param_name, value in scalar_params.items():
            val = value.value if isinstance(value, RefValue) else value
            df_state.loc[grid_id, (param_name, "0")] = val

        df_hourly_profile = self.snowprof_24hr.to_df_state(grid_id, "snowprof_24hr")
        df_state = df_state.combine_first(df_hourly_profile)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "SnowParams":
        """
        Reconstruct SnowParams from a DataFrame state format.

        Args:
            df: DataFrame containing snow parameters.
            grid_id: Grid ID for the DataFrame index.

        Returns:
            SnowParams: Instance of SnowParams.
        """
        # Extract scalar attributes
        scalar_params = {
            "crwmax": df.loc[grid_id, ("crwmax", "0")],
            "crwmin": df.loc[grid_id, ("crwmin", "0")],
            "narp_emis_snow": df.loc[grid_id, ("narp_emis_snow", "0")],
            "preciplimit": df.loc[grid_id, ("preciplimit", "0")],
            "preciplimitalb": df.loc[grid_id, ("preciplimitalb", "0")],
            "snowalbmax": df.loc[grid_id, ("snowalbmax", "0")],
            "snowalbmin": df.loc[grid_id, ("snowalbmin", "0")],
            "snowdensmin": df.loc[grid_id, ("snowdensmin", "0")],
            "snowdensmax": df.loc[grid_id, ("snowdensmax", "0")],
            "snowlimbldg": df.loc[grid_id, ("snowlimbldg", "0")],
            "snowlimpaved": df.loc[grid_id, ("snowlimpaved", "0")],
            "tau_a": df.loc[grid_id, ("tau_a", "0")],
            "tau_f": df.loc[grid_id, ("tau_f", "0")],
            "tau_r": df.loc[grid_id, ("tau_r", "0")],
            "tempmeltfact": df.loc[grid_id, ("tempmeltfact", "0")],
            "radmeltfact": df.loc[grid_id, ("radmeltfact", "0")],
        }

        # Convert scalar parameters to RefValue
        scalar_params = {
            key: RefValue(value) for key, value in scalar_params.items()
        }

        # Extract HourlyProfile
        snowprof_24hr = HourlyProfile.from_df_state(df, grid_id, "snowprof_24hr")

        # Construct and return the SnowParams instance
        return cls(snowprof_24hr=snowprof_24hr, **scalar_params)


class LandCover(BaseModel):
    paved: PavedProperties = Field(
        default_factory=PavedProperties,
        description="Properties for paved surfaces like roads and pavements",
    )
    bldgs: BldgsProperties = Field(
        default_factory=BldgsProperties,
        description="Properties for building surfaces including roofs and walls",
    )
    evetr: EvetrProperties = Field(
        default_factory=EvetrProperties,
        description="Properties for evergreen trees and vegetation",
    )
    dectr: DectrProperties = Field(
        default_factory=DectrProperties,
        description="Properties for deciduous trees and vegetation",
    )
    grass: GrassProperties = Field(
        default_factory=GrassProperties, description="Properties for grass surfaces"
    )
    bsoil: BsoilProperties = Field(
        default_factory=BsoilProperties, description="Properties for bare soil surfaces"
    )
    water: WaterProperties = Field(
        default_factory=WaterProperties,
        description="Properties for water surfaces like lakes and ponds",
    )

    ref: Optional[Reference] = None

    @model_validator(mode="after")
    def set_surface_types(self) -> "LandCover":
        # Set surface types and validate
        surface_map = {
            "paved": (self.paved, SurfaceType.PAVED),
            "bldgs": (self.bldgs, SurfaceType.BLDGS),
            "dectr": (self.dectr, SurfaceType.DECTR),
            "evetr": (self.evetr, SurfaceType.EVETR),
            "grass": (self.grass, SurfaceType.GRASS),
            "bsoil": (self.bsoil, SurfaceType.BSOIL),
            "water": (self.water, SurfaceType.WATER),
        }

        for prop, surface_type in surface_map.values():
            prop.set_surface_type(surface_type)

        return self

    @model_validator(mode="after")
    def validate_land_cover_fractions(self) -> "LandCover":
        fractions = {
            "paved": self.paved.sfr.value,
            "bldgs": self.bldgs.sfr.value,
            "evetr": self.evetr.sfr.value,
            "dectr": self.dectr.sfr.value,
            "grass": self.grass.sfr.value,
            "bsoil": self.bsoil.sfr.value,
            "water": self.water.sfr.value,
        }

        total = sum(fractions.values())
        if abs(total - 1.0) > 1e-6:
            details = ", ".join(f"{k}={v:.3f}" for k, v in fractions.items())
            raise ValueError(f"Land cover fractions must sum to 1.0 (got {total:.6f}): {details}")
        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert land cover to DataFrame state format"""
        # df_state = init_df_state(grid_id)

        list_df_state = []
        for lc in ["paved", "bldgs", "dectr", "evetr", "grass", "bsoil", "water"]:
            df_state = getattr(self, lc).to_df_state(grid_id)
            list_df_state.append(df_state)
        df_state = pd.concat(list_df_state, axis=1)
        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "LandCover":
        """Reconstruct LandCover instance from DataFrame state.

        Args:
            df: DataFrame containing land cover parameters
            grid_id: Grid ID for the DataFrame index

        Returns:
            LandCover: Reconstructed LandCover instance
        """
        # Reconstruct each surface type from the DataFrame
        params = {
            "paved": PavedProperties.from_df_state(df, grid_id),
            "bldgs": BldgsProperties.from_df_state(df, grid_id),
            "evetr": EvetrProperties.from_df_state(df, grid_id),
            "dectr": DectrProperties.from_df_state(df, grid_id),
            "grass": GrassProperties.from_df_state(df, grid_id),
            "bsoil": BsoilProperties.from_df_state(df, grid_id),
            "water": WaterProperties.from_df_state(df, grid_id),
        }

        # Return reconstructed instance
        return cls(**params)


class ArchetypeProperties(BaseModel):
    # Not used in STEBBS - DAVE only
    # BuildingCode='1'
    # BuildingClass='SampleClass'

    BuildingType: str = "SampleType"
    BuildingName: str = "SampleBuilding"
    BuildingCount: FlexibleRefValue(int) = Field(
        default=1, description="Number of buildings of this archetype [-]", json_schema_extra={"unit": "dimensionless"}
    )
    Occupants: FlexibleRefValue(int) = Field(
        default=1,
        description="Number of occupants present in building [-]", json_schema_extra={"unit": "dimensionless"},
    )

    # Not used in STEBBS - DAVE only
    # hhs0: int = Field(default=0, description="")
    # hhs1: int = Field(default=0, description="")
    # hhs2: int = Field(default=0, description="")
    # hhs3: int = Field(default=0, description="")
    # hhs4: int = Field(default=0, description="")
    # hhs5: int = Field(default=0, description="")
    # hhs6: int = Field(default=0, description="")
    # hhs7: int = Field(default=0, description="")
    # hhs8: int = Field(default=0, description="")
    # age_0_4: int = Field(default=0, description="")
    # age_5_11: int = Field(default=0, description="")
    # age_12_18: int = Field(default=0, description="")
    # age_19_64: int = Field(default=0, description="")
    # age_65plus: int = Field(default=0, description="")

    stebbs_Height: FlexibleRefValue(float) = Field(
        default=10.0,
        description="Building height [m]", json_schema_extra={"unit": "m"},
        gt=0.0,
    )
    FootprintArea: FlexibleRefValue(float) = Field(
        default=64.0,
        description="Building footprint area [m2]", json_schema_extra={"unit": "m^2"},
        gt=0.0,
    )
    WallExternalArea: FlexibleRefValue(float) = Field(
        default=80.0,
        description="External wall area (including window area) [m2]",
        json_schema_extra={"unit": "m^2"},
        gt=0.0,
    )
    RatioInternalVolume: FlexibleRefValue(float) = Field(
        default=0.01,
        description="Ratio of internal mass volume to total building volume [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WWR: FlexibleRefValue(float) = Field(
        default=0.20,
        description="window to wall ratio [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WallThickness: FlexibleRefValue(float) = Field(
        default=20.0,
        description="Thickness of external wall and roof (weighted) [m]",
        json_schema_extra={"unit": "m"},
        gt=0.0,
    )
    WallEffectiveConductivity: FlexibleRefValue(float) = Field(
        default=60.0,
        description="Effective thermal conductivity of walls and roofs (weighted) [W m-1 K-1]",
        json_schema_extra={"unit": "W m^-1 K^-1"},
        gt=0.0,
    )
    WallDensity: FlexibleRefValue(float) = Field(
        default=1600.0,
        description="Effective density of the walls and roof (weighted) [kg m-3]",
        json_schema_extra={"unit": "kg m^-3"},
        gt=0.0,
    )
    WallCp: FlexibleRefValue(float) = Field(
        default=850.0,
        description="Effective specific heat capacity of walls and roof (weighted) [J kg-1 K-1]",
        json_schema_extra={"unit": "J kg^-1 K^-1"},
        gt=0.0,
    )
    Wallx1: FlexibleRefValue(float) = Field(
        default=1.0,
        description="Weighting factor for heat capacity of walls and roof [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WallExternalEmissivity: FlexibleRefValue(float) = Field(
        default=0.9,
        description="Emissivity of the external surface of walls and roof [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WallInternalEmissivity: FlexibleRefValue(float) = Field(
        default=0.9,
        description="Emissivity of the internal surface of walls and roof [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WallTransmissivity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Transmissivity of walls and roof [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WallAbsorbtivity: FlexibleRefValue(float) = Field(
        default=0.8,
        description="Absorbtivity of walls and roof [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WallReflectivity: FlexibleRefValue(float) = Field(
        default=0.2,
        description="Reflectivity of the external surface of walls and roof [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    FloorThickness: FlexibleRefValue(float) = Field(
        default=0.2,
        description="Thickness of ground floor [m]", json_schema_extra={"unit": "m"},
        gt=0.0,
    )
    GroundFloorEffectiveConductivity: FlexibleRefValue(float) = Field(
        default=0.15,
        description="Effective thermal conductivity of ground floor [W m-1 K-1]", json_schema_extra={"unit": "W m^-1 K^-1"},
        gt=0.0,
    )
    GroundFloorDensity: FlexibleRefValue(float) = Field(
        default=500.0,
        description="Density of the ground floor [kg m-3]", json_schema_extra={"unit": "kg m^-3"},
        gt=0.0,
    )
    GroundFloorCp: FlexibleRefValue(float) = Field(
        default=1500.0,
        description="Effective specific heat capacity of the ground floor [J kg-1 K-1]", json_schema_extra={"unit": "J kg^-1 K^-1"},
        gt=0.0,
    )
    WindowThickness: FlexibleRefValue(float) = Field(
        default=0.015,
        description="Window thickness [m]", json_schema_extra={"unit": "m"},
        gt=0.0,
    )
    WindowEffectiveConductivity: FlexibleRefValue(float) = Field(
        default=1.0,
        description="Effective thermal conductivity of windows [W m-1 K-1]", json_schema_extra={"unit": "W m^-1 K^-1"},
        gt=0.0,
    )
    WindowDensity: FlexibleRefValue(float) = Field(
        default=2500.0,
        description="Effective density of the windows [kg m-3]", json_schema_extra={"unit": "kg m^-3"},
        gt=0.0,
    )
    WindowCp: FlexibleRefValue(float) = Field(
        default=840.0,
        description="Effective specific heat capacity of windows [J kg-1 K-1]", json_schema_extra={"unit": "J kg^-1 K^-1"},
        gt=0.0,
    )
    WindowExternalEmissivity: FlexibleRefValue(float) = Field(
        default=0.90,
        description="Emissivity of the external surface of windows [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WindowInternalEmissivity: FlexibleRefValue(float) = Field(
        default=0.90,
        description="Emissivity of the internal surface of windows [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WindowTransmissivity: FlexibleRefValue(float) = Field(
        default=0.90,
        description="Transmissivity of windows [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WindowAbsorbtivity: FlexibleRefValue(float) = Field(
        default=0.01,
        description="Absorbtivity of windows [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    WindowReflectivity: FlexibleRefValue(float) = Field(
        default=0.09,
        description="Reflectivity of the external surface of windows [-]", json_schema_extra={"unit": "dimensionless"},
        ge=0.0,
        le=1.0,
    )
    # TODO: Add defaults below here
    InternalMassDensity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Effective density of the internal mass [kg m-3]", json_schema_extra={"unit": "kg m^-3"},
    )
    InternalMassCp: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Specific heat capacity of internal mass [J kg-1 K-1]", json_schema_extra={"unit": "J kg^-1 K^-1"},
    )
    InternalMassEmissivity: FlexibleRefValue(float) = Field(
        default=0.0, description="Emissivity of internal mass [-]", json_schema_extra={"unit": "dimensionless"},
    )
    MaxHeatingPower: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Maximum power demand of heating system [W]", json_schema_extra={"unit": "W"},
    )
    WaterTankWaterVolume: FlexibleRefValue(float) = Field(
        default=0.0, description="Volume of water in hot water tank [m3]", json_schema_extra={"unit": "m^3"},
    )
    MaximumHotWaterHeatingPower: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Maximum power demand of water heating system [W]", json_schema_extra={"unit": "W"},
    )
    HeatingSetpointTemperature: FlexibleRefValue(float) = Field(
        default=0.0, description="Heating setpoint temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    CoolingSetpointTemperature: FlexibleRefValue(float) = Field(
        default=0.0, description="Cooling setpoint temperature [degC]", json_schema_extra={"unit": "degC"},
    )

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert ArchetypeProperties to DataFrame state format."""

        df_state = init_df_state(grid_id)

        # Create an empty DataFrame with MultiIndex columns
        columns = [
            (field.lower(), "0") for field in self.__class__.model_fields.keys() if field != "ref"
        ]
        df_state = pd.DataFrame(
            index=[grid_id], columns=pd.MultiIndex.from_tuples(columns)
        )

        # Set the values in the DataFrame
        for field_name, field_info in self.__class__.model_fields.items():
            if field_name == "ref":
                continue
            attribute = getattr(self, field_name)
            if isinstance(attribute, RefValue):
                value = attribute.value
            else:
                value = attribute
            df_state.loc[grid_id, (field_name.lower(), "0")] = value

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "ArchetypeProperties":
        """Reconstruct ArchetypeProperties from DataFrame state format."""
        # Extract the values from the DataFrame
        params = {
            field_name: df.loc[grid_id, (field_name.lower(), "0")]
            for field_name in cls.model_fields.keys()
            if field_name != "ref"
        }

        # Convert params to RefValue
        non_value_with_doi = ["BuildingType", "BuildingName"]
        params = {
            key: (RefValue(value) if key not in non_value_with_doi else value)
            for key, value in params.items()
        }

        # Create an instance using the extracted parameters
        return cls(**params)


class StebbsProperties(BaseModel):
    WallInternalConvectionCoefficient: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Internal convection coefficient of walls and roof [W m-2 K-1]", json_schema_extra={"unit": "W m^-2 K^-1"},
    )
    InternalMassConvectionCoefficient: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Convection coefficient of internal mass [W m-2 K-1]", json_schema_extra={"unit": "W m^-2 K^-1"},
    )
    FloorInternalConvectionCoefficient: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Internal convection coefficient of ground floor [W m-2 K-1]", json_schema_extra={"unit": "W m^-2 K^-1"},
    )
    WindowInternalConvectionCoefficient: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Internal convection coefficient of windows [W m-2 K-1]", json_schema_extra={"unit": "W m^-2 K^-1"},
    )
    WallExternalConvectionCoefficient: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial external convection coefficient of walls and roof [W m-2 K-1]", json_schema_extra={"unit": "W m^-2 K^-1"},
    )
    WindowExternalConvectionCoefficient: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial external convection coefficient of windows [W m-2 K-1]", json_schema_extra={"unit": "W m^-2 K^-1"},
    )
    GroundDepth: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Depth of external ground (deep soil) [m]",
        json_schema_extra={"unit": "m"}
    )
    ExternalGroundConductivity: FlexibleRefValue(float) = Field(
        default=0.0, description="External ground thermal conductivity", json_schema_extra={"unit": "W m^-1 K^-1"}
    )
    IndoorAirDensity: FlexibleRefValue(float) = Field(
        default=0.0, description="Density of indoor air [kg m-3]", json_schema_extra={"unit": "kg m^-3"}
    )
    IndoorAirCp: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Specific heat capacity of indoor air [J kg-1 K-1]", json_schema_extra={"unit": "J kg^-1 K^-1"},
    )
    WallBuildingViewFactor: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Building view factor of external walls [-]", json_schema_extra={"unit": "dimensionless"},
    )
    WallGroundViewFactor: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Ground view factor of external walls [-]", json_schema_extra={"unit": "dimensionless"},
    )
    WallSkyViewFactor: FlexibleRefValue(float) = Field(
        default=0.0, description="Sky view factor of external walls [-]", json_schema_extra={"unit": "dimensionless"}
    )
    MetabolicRate: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Metabolic rate of building occupants [W]", json_schema_extra={"unit": "W"},
    )
    LatentSensibleRatio: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Latent-to-sensible ratio of metabolic energy release of occupants [-]", json_schema_extra={"unit": "dimensionless"},
    )
    ApplianceRating: FlexibleRefValue(float) = Field(
        default=0.0, description="Power demand of single appliance [W]", json_schema_extra={"unit": "W"}
    )
    TotalNumberofAppliances: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Number of appliances present in building [-]", json_schema_extra={"unit": "dimensionless"},
    )
    ApplianceUsageFactor: FlexibleRefValue(float) = Field(
        default=0.0, description="Number of appliances in use [-]", json_schema_extra={"unit": "dimensionless"}
    )
    HeatingSystemEfficiency: FlexibleRefValue(float) = Field(
        default=0.0, description="Efficiency of space heating system [-]", json_schema_extra={"unit": "dimensionless"}
    )
    MaxCoolingPower: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Maximum power demand of cooling system [W]", json_schema_extra={"unit": "W"},
    )
    CoolingSystemCOP: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Coefficient of performance of cooling system [-]", json_schema_extra={"unit": "dimensionless"},
    )
    VentilationRate: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Ventilation rate (air changes per hour, ACH) [h-1]",
        json_schema_extra={"unit": "h^-1"}
    )
    IndoorAirStartTemperature: FlexibleRefValue(float) = Field(
        default=0.0, description="Initial indoor air temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    IndoorMassStartTemperature: FlexibleRefValue(float) = Field(
        default=0.0, description="Initial indoor mass temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    WallIndoorSurfaceTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial wall/roof indoor surface temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    WallOutdoorSurfaceTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial wall/roof outdoor surface temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    WindowIndoorSurfaceTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial window indoor surface temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    WindowOutdoorSurfaceTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial window outdoor surface temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    GroundFloorIndoorSurfaceTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial ground floor indoor surface temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    GroundFloorOutdoorSurfaceTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial ground floor outdoor surface temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    WaterTankTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial water temperature in hot water tank [degC]", json_schema_extra={"unit": "degC"},
    )
    InternalWallWaterTankTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial hot water tank internal wall temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    ExternalWallWaterTankTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial hot water tank external wall temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    WaterTankWallThickness: FlexibleRefValue(float) = Field(
        default=0.0, description="Hot water tank wall thickness [m]", json_schema_extra={"unit": "m"},
    )
    MainsWaterTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Temperature of water coming into the water tank [degC]", json_schema_extra={"unit": "degC"},
    )
    WaterTankSurfaceArea: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Surface area of hot water tank cylinder [m2]", json_schema_extra={"unit": "m^2"},
    )
    HotWaterHeatingSetpointTemperature: FlexibleRefValue(float) = Field(
        default=0.0, description="Water tank setpoint temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    HotWaterTankWallEmissivity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Effective external wall emissivity of the hot water tank [-]", json_schema_extra={"unit": "dimensionless"},
    )
    DomesticHotWaterTemperatureInUseInBuilding: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial water temperature of water held in use in building [degC]", json_schema_extra={"unit": "degC"},
    )
    InternalWallDHWVesselTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial hot water vessel internal wall temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    ExternalWallDHWVesselTemperature: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Initial hot water vessel external wall temperature [degC]", json_schema_extra={"unit": "degC"},
    )
    DHWVesselWallThickness: FlexibleRefValue(float) = Field(
        default=0.0, description="Hot water vessel wall thickness [m]", json_schema_extra={"unit": "m"},
    )
    DHWWaterVolume: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Volume of water held in use in building [m3]", json_schema_extra={"unit": "m^3"},
    )
    DHWSurfaceArea: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Surface area of hot water in vessels in building [m2]", json_schema_extra={"unit": "m^2"},
    )
    DHWVesselEmissivity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="NEEDS CHECKED! NOT USED (assumed same as DHWVesselWallEmissivity) [-]",
        json_schema_extra={"unit": "dimensionless"}
    )
    HotWaterFlowRate: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Hot water flow rate from tank to vessel [m3 s-1]", json_schema_extra={"unit": "m^3 s^-1"},
    )
    DHWDrainFlowRate: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Flow rate of hot water held in building to drain [m3 s-1]", json_schema_extra={"unit": "m^3 s^-1"},
    )
    DHWSpecificHeatCapacity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Specific heat capacity of hot water [J kg-1 K-1]", json_schema_extra={"unit": "J kg^-1 K^-1"},
    )
    HotWaterTankSpecificHeatCapacity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Specific heat capacity of hot water tank wal [J kg-1 K-1]", json_schema_extra={"unit": "J kg^-1 K^-1"},
    )
    DHWVesselSpecificHeatCapacity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Specific heat capacity of vessels containing hot water in use in buildings [J kg-1 K-1]", json_schema_extra={"unit": "J kg^-1 K^-1"},
    )
    DHWDensity: FlexibleRefValue(float) = Field(
        default=0.0, description="Density of hot water in use [kg m-3]", json_schema_extra={"unit": "kg m^-3"},
    )
    HotWaterTankWallDensity: FlexibleRefValue(float) = Field(
        default=0.0, description="Density of hot water tank wall [kg m-3]", json_schema_extra={"unit": "kg m^-3"},
    )
    DHWVesselDensity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Density of vessels containing hot water in use [kg m-3]", json_schema_extra={"unit": "kg m^-3"},
    )
    HotWaterTankBuildingWallViewFactor: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Water tank/vessel internal building wall/roof view factor [-]", json_schema_extra={"unit": "dimensionless"},
    )
    HotWaterTankInternalMassViewFactor: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Water tank/vessel building internal mass view factor [-]", json_schema_extra={"unit": "dimensionless"},
    )
    HotWaterTankWallConductivity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Effective wall conductivity of the hot water tank [W m-1 K-1]", json_schema_extra={"unit": "W m^-1 K^-1"},
    )
    HotWaterTankInternalWallConvectionCoefficient: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Effective internal wall convection coefficient of the hot water tank [W m-2 K-1]", json_schema_extra={"unit": "W m^-2 K^-1"},
    )
    HotWaterTankExternalWallConvectionCoefficient: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Effective external wall convection coefficient of the hot water tank [W m-2 K-1]", json_schema_extra={"unit": "W m^-2 K^-1"},
    )
    DHWVesselWallConductivity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Effective wall conductivity of the hot water tank [W m-1 K-1]", json_schema_extra={"unit": "W m^-1 K^-1"},
    )
    DHWVesselInternalWallConvectionCoefficient: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Effective internal wall convection coefficient of the vessels holding hot water in use in building [W m-2 K-1]", json_schema_extra={"unit": "W m^-2 K^-1"},
    )
    DHWVesselExternalWallConvectionCoefficient: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Effective external wall convection coefficient of the vessels holding hot water in use in building [W m-2 K-1]", json_schema_extra={"unit": "W m^-2 K^-1"},
    )
    DHWVesselWallEmissivity: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Effective external wall emissivity of hot water being used within building [-]", json_schema_extra={"unit": "dimensionless"},
    )
    HotWaterHeatingEfficiency: FlexibleRefValue(float) = Field(
        default=0.0, description="Efficiency of hot water system [-]", json_schema_extra={"unit": "dimensionless"}
    )
    MinimumVolumeOfDHWinUse: FlexibleRefValue(float) = Field(
        default=0.0, description="Minimum volume of hot water in use [m3]", json_schema_extra={"unit": "m^3"}
    )

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert StebbsProperties to DataFrame state format."""
        df_state = init_df_state(grid_id)

        # Create an empty DataFrame with MultiIndex columns
        columns = [
            (field.lower(), "0") for field in self.__class__.model_fields.keys() if field != "ref"
        ]
        df_state = pd.DataFrame(
            index=[grid_id], columns=pd.MultiIndex.from_tuples(columns)
        )

        # Set the values in the DataFrame
        for field_name, field_info in self.__class__.model_fields.items():
            if field_name == "ref":
                continue
            field_val = getattr(self, field_name)
            val = field_val.value if isinstance(field_val, RefValue) else field_val
            df_state.loc[grid_id, (field_name.lower(), "0")] = val

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "StebbsProperties":
        """Reconstruct StebbsProperties from DataFrame state format."""
        # Extract the values from the DataFrame
        params = {
            field_name: df.loc[grid_id, (field_name.lower(), "0")]
            for field_name in cls.model_fields.keys()
            if field_name != "ref"
        }

        # Convert params to RefValue
        params = {key: RefValue(value) for key, value in params.items()}

        # Create an instance using the extracted parameters
        return cls(**params)


class SPARTACUSParams(BaseModel):
    air_ext_lw: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Air extinction coefficient for longwave radiation", json_schema_extra={"unit": "m^-1"},
    )
    air_ext_sw: FlexibleRefValue(float) = Field(
        default=0.0,
        description="Air extinction coefficient for shortwave radiation", json_schema_extra={"unit": "m^-1"},
    )
    air_ssa_lw: FlexibleRefValue(float) = Field(
        default=0.5,
        description="Air single scattering albedo for longwave radiation", json_schema_extra={"unit": "dimensionless"},
    )
    air_ssa_sw: FlexibleRefValue(float) = Field(
        default=0.5,
        description="Air single scattering albedo for shortwave radiation", json_schema_extra={"unit": "dimensionless"},
    )
    ground_albedo_dir_mult_fact: FlexibleRefValue(float) = Field(
        default=1.0,
        description="Multiplication factor for direct ground albedo", json_schema_extra={"unit": "dimensionless"}
    )
    n_stream_lw_urban: FlexibleRefValue(int) = Field(
        default=2,
        description="Number of streams for longwave radiation in urban areas", json_schema_extra={"unit": "dimensionless"},
    )
    n_stream_sw_urban: FlexibleRefValue(int) = Field(
        default=2,
        description="Number of streams for shortwave radiation in urban areas", json_schema_extra={"unit": "dimensionless"},
    )
    n_vegetation_region_urban: FlexibleRefValue(int) = Field(
        default=1,
        description="Number of vegetation regions in urban areas", json_schema_extra={"unit": "dimensionless"},
    )
    sw_dn_direct_frac: FlexibleRefValue(float) = Field(
        default=0.5,
        description="Fraction of downward shortwave radiation that is direct", json_schema_extra={"unit": "dimensionless"},
    )
    use_sw_direct_albedo: FlexibleRefValue(float) = Field(
        default=1.0,
        description="Flag to use direct albedo for shortwave radiation", json_schema_extra={"unit": "dimensionless"}
    )
    veg_contact_fraction_const: FlexibleRefValue(float) = Field(
        default=0.5, description="Constant vegetation contact fraction", json_schema_extra={"unit": "dimensionless"}
    )
    veg_fsd_const: FlexibleRefValue(float) = Field(
        default=0.5,
        description="Constant vegetation fractional standard deviation", json_schema_extra={"unit": "dimensionless"},
    )
    veg_ssa_lw: FlexibleRefValue(float) = Field(
        default=0.5,
        description="Vegetation single scattering albedo for longwave radiation", json_schema_extra={"unit": "dimensionless"},
    )
    veg_ssa_sw: FlexibleRefValue(float) = Field(
        default=0.5,
        description="Vegetation single scattering albedo for shortwave radiation", json_schema_extra={"unit": "dimensionless"},
    )

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """
        Convert SPARTACUS parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing SPARTACUS parameters
        """
        # Initialize DataFrame with grid index
        df_state = init_df_state(grid_id)

        # Map SPARTACUS parameters to DataFrame columns
        spartacus_params = {
            "air_ext_lw": self.air_ext_lw,
            "air_ext_sw": self.air_ext_sw,
            "air_ssa_lw": self.air_ssa_lw,
            "air_ssa_sw": self.air_ssa_sw,
            "ground_albedo_dir_mult_fact": self.ground_albedo_dir_mult_fact,
            "n_stream_lw_urban": self.n_stream_lw_urban,
            "n_stream_sw_urban": self.n_stream_sw_urban,
            "n_vegetation_region_urban": self.n_vegetation_region_urban,
            "sw_dn_direct_frac": self.sw_dn_direct_frac,
            "use_sw_direct_albedo": self.use_sw_direct_albedo,
            "veg_contact_fraction_const": self.veg_contact_fraction_const,
            "veg_fsd_const": self.veg_fsd_const,
            "veg_ssa_lw": self.veg_ssa_lw,
            "veg_ssa_sw": self.veg_ssa_sw,
        }

        # Assign each parameter to its corresponding column in the DataFrame
        for param_name, value in spartacus_params.items():
            val = value.value if isinstance(value, RefValue) else value
            df_state[(param_name, "0")] = val

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "SPARTACUSParams":
        """
        Reconstruct SPARTACUSParams from DataFrame state format.

        Args:
            df: DataFrame containing SPARTACUS parameters
            grid_id: Grid ID for the DataFrame index

        Returns:
            SPARTACUSParams: An instance of SPARTACUSParams
        """

        spartacus_params = {
            "air_ext_lw",
            "air_ext_sw",
            "air_ssa_lw",
            "air_ssa_sw",
            "ground_albedo_dir_mult_fact",
            "n_stream_lw_urban",
            "n_stream_sw_urban",
            "n_vegetation_region_urban",
            "sw_dn_direct_frac",
            "use_sw_direct_albedo",
            "veg_contact_fraction_const",
            "veg_fsd_const",
            "veg_ssa_lw",
            "veg_ssa_sw",
        }

        params = {
            param: RefValue(df.loc[grid_id, (param, "0")])
            for param in spartacus_params
        }

        return cls(**params)


class LUMPSParams(BaseModel):
    raincover: FlexibleRefValue(float) = Field(
        ge=0, le=1, default=0.25, 
        description="Rain water coverage fraction", 
        json_schema_extra={"unit": "dimensionless"}
    )
    rainmaxres: FlexibleRefValue(float) = Field(
        ge=0, le=20, default=0.25, 
        description="Maximum rain water storage", 
        json_schema_extra={"unit": "mm"}
    )
    drainrt: FlexibleRefValue(float) = Field(
        ge=0, le=1, default=0.25, 
        description="Drainage rate coefficient", 
        json_schema_extra={"unit": "dimensionless"}
    )
    veg_type: FlexibleRefValue(int) = Field(
        default=1, 
        description="Vegetation type selection", 
        json_schema_extra={"unit": "dimensionless"}
    )

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert LUMPS parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing LUMPS parameters
        """
        df_state = init_df_state(grid_id)

        # Add all attributes
        for attr in ["raincover", "rainmaxres", "drainrt", "veg_type"]:
            field_val = getattr(self, attr)
            val = field_val.value if isinstance(field_val, RefValue) else field_val
            df_state[(attr, "0")] = val

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "LUMPSParams":
        """Create LUMPSParams from DataFrame state format.

        Args:
            df: DataFrame containing LUMPS parameters
            grid_id: Grid ID for the DataFrame index

        Returns:
            LUMPSParams: Instance of LUMPSParams
        """
        # Extract attributes from DataFrame
        params = {}
        for attr in ["raincover", "rainmaxres", "drainrt", "veg_type"]:
            params[attr] = df.loc[grid_id, (attr, "0")]

        # Convert attributes to RefValue
        params = {key: RefValue(value) for key, value in params.items()}

        return cls(**params)


class SiteProperties(BaseModel):
    lat: FlexibleRefValue(float) = Field(
        ge=-90,
        le=90,
        description="Latitude of the site in degrees", json_schema_extra={"unit": "degrees"},
        default=51.5,
    )
    lng: FlexibleRefValue(float) = Field(
        ge=-180,
        le=180,
        description="Longitude of the site in degrees", json_schema_extra={"unit": "degrees"},
        default=-0.13,
    )
    alt: FlexibleRefValue(float) = Field(
        gt=0,
        description="Altitude of the site above sea level", json_schema_extra={"unit": "m"},
        default=40.0,
    )
    timezone: FlexibleRefValue(int) = Field(
        ge=-12,
        le=12,
        description="Time zone offset from UTC", json_schema_extra={"unit": "hour"},
        default=0,
    )
    surfacearea: FlexibleRefValue(float) = Field(
        gt=0,
        description="Total surface area of the site",
        json_schema_extra={"unit": "m^2"},
        default=10000.0,  # 1 hectare in m²
    )
    z: FlexibleRefValue(float) = Field(
        gt=0,
        description="Measurement height", json_schema_extra={"unit": "m"},
        default=10.0
    )
    z0m_in: FlexibleRefValue(float) = Field(
        gt=0,
        description="Momentum roughness length", json_schema_extra={"unit": "m"},
        default=1.0,
    )
    zdm_in: FlexibleRefValue(float) = Field(
        gt=0,
        description="Zero-plane displacement height", json_schema_extra={"unit": "m"},
        default=5.0,
    )
    pipecapacity: FlexibleRefValue(float) = Field(
        gt=0,
        description="Maximum capacity of drainage pipes", json_schema_extra={"unit": "mm h^-1"},
        default=100.0,
    )
    runofftowater: FlexibleRefValue(float) = Field(
        ge=0,
        le=1,
        description="Fraction of excess water going to water bodies", json_schema_extra={"unit": "dimensionless"},
        default=0.0,
    )
    narp_trans_site: FlexibleRefValue(float) = Field(
        description="Site-specific NARP transmission coefficient", json_schema_extra={"unit": "dimensionless"},
        default=0.2,
    )
    lumps: LUMPSParams = Field(
        default_factory=LUMPSParams,
        description="Parameters for Local-scale Urban Meteorological Parameterization Scheme",
    )
    spartacus: SPARTACUSParams = Field(
        default_factory=SPARTACUSParams,
        description="Parameters for Solar Parametrizations for Radiative Transfer through Urban Canopy Scheme",
    )
    stebbs: StebbsProperties = Field(
        default_factory=StebbsProperties,
        description="Parameters for the STEBBS building energy model",
    )
    building_archetype: ArchetypeProperties = Field(
        default_factory=ArchetypeProperties,
        description="Parameters for building archetypes",
    )
    conductance: Conductance = Field(
        default_factory=Conductance,
        description="Parameters for surface conductance calculations",
    )
    irrigation: IrrigationParams = Field(
        default_factory=IrrigationParams,
        description="Parameters for irrigation modelling",
    )
    anthropogenic_emissions: AnthropogenicEmissions = Field(
        default_factory=AnthropogenicEmissions,
        description="Parameters for anthropogenic heat and water emissions",
    )
    snow: SnowParams = Field(
        default_factory=SnowParams, description="Parameters for snow modelling"
    )
    land_cover: LandCover = Field(
        default_factory=LandCover,
        description="Parameters for land cover characteristics",
    )
    vertical_layers: VerticalLayers = Field(
        default_factory=VerticalLayers,
        description="Parameters for vertical layer structure",
    )

    n_buildings: FlexibleRefValue(int) = Field(
        default=1,
        description="Number of buildings in the site", json_schema_extra={"unit": "dimensionless"},
    )

    h_std: FlexibleRefValue(float) = Field(
        default=10.0,
        description="Standard deviation of building heights in the site", json_schema_extra={"unit": "m"},
    )

    lambda_c: FlexibleRefValue(float) = Field(
        default=0,
        description="External building surface area to plan area ratio", json_schema_extra={"unit": "m^2 m^-2"},
        ge=0
    )

    ref: Optional[Reference] = None

    model_config = ConfigDict(
        extra="forbid",  # This will prevent extra fields from being accepted
        validate_assignment=True,  # This will validate fields on assignment
        validate_default=True,  # This will validate default values
    )

    @model_validator(mode="after")
    def validate_required_fields(self) -> "SiteProperties":
        """Validate that all required fields are present and have valid values."""
        errors = []

        # List of required fields that must be present and non-None
        required_fields = [
            "lat",
            "lng",
            "alt",
            "timezone",
            "surfacearea",
            "z",
            "z0m_in",
            "zdm_in",
            "pipecapacity",
            "runofftowater",
            "narp_trans_site",
            "lumps",
            "spartacus",
            "conductance",
            "irrigation",
            "anthropogenic_emissions",
            "snow",
            "land_cover",
            "vertical_layers",
        ]

        for field in required_fields:
            value = getattr(self, field, None)
            if value is None:
                errors.append(f"Required field '{field}' is missing")
            elif isinstance(value, RefValue) and (value.value if isinstance(value, RefValue) else value) is None:
                errors.append(f"Required field '{field}' has no value")

        # Additional validation rules
        z0m_val = self.z0m_in.value if isinstance(self.z0m_in, RefValue) else self.z0m_in
        zdm_val = self.zdm_in.value if isinstance(self.zdm_in, RefValue) else self.zdm_in
        if z0m_val >= zdm_val:
            errors.append(
                f"z0m_in ({z0m_val}) must be less than zdm_in ({zdm_val})"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert site properties to DataFrame state format"""
        df_state = init_df_state(grid_id)

        # simple attributes
        for var in [
            "lat",
            "lng",
            "alt",
            "timezone",
            "surfacearea",
            "z",
            "z0m_in",
            "zdm_in",
            "pipecapacity",
            "runofftowater",
            "narp_trans_site",
            "n_buildings",
            "h_std",
            "lambda_c"
        ]:
            field_val = getattr(self, var)
            val = field_val.value if isinstance(field_val, RefValue) else field_val
            df_state.loc[grid_id, (f"{var}", "0")] = val

        # complex attributes
        df_lumps = self.lumps.to_df_state(grid_id)
        df_spartacus = self.spartacus.to_df_state(grid_id)
        df_conductance = self.conductance.to_df_state(grid_id)
        df_irrigation = self.irrigation.to_df_state(grid_id)
        df_anthropogenic_emissions = self.anthropogenic_emissions.to_df_state(grid_id)
        df_snow = self.snow.to_df_state(grid_id)
        df_land_cover = self.land_cover.to_df_state(grid_id)
        df_vertical_layers = self.vertical_layers.to_df_state(grid_id)
        df_stebbs = self.stebbs.to_df_state(grid_id)
        df_building_archetype = self.building_archetype.to_df_state(grid_id)

        df_state = pd.concat(
            [
                df_state,
                df_lumps,
                df_spartacus,
                df_conductance,
                df_irrigation,
                df_anthropogenic_emissions,
                df_snow,
                df_land_cover,
                df_vertical_layers,
                df_stebbs,
                df_building_archetype,
            ],
            axis=1,
        )
        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "SiteProperties":
        """Reconstruct SiteProperties from DataFrame state format.

        Args:
            df: DataFrame containing site properties
            grid_id: Grid ID for the DataFrame index

        Returns:
            SiteProperties: Reconstructed instance
        """
        # Extract simple attributes
        params = {}
        for var in [
            "lat",
            "lng",
            "alt",
            "timezone",
            "surfacearea",
            "z",
            "z0m_in",
            "zdm_in",
            "pipecapacity",
            "runofftowater",
            "narp_trans_site",
            "n_buildings",
            "h_std",
            "lambda_c"
        ]:
            params[var] = RefValue(df.loc[grid_id, (var, "0")])

        # Extract complex attributes
        params["lumps"] = LUMPSParams.from_df_state(df, grid_id)
        params["spartacus"] = SPARTACUSParams.from_df_state(df, grid_id)
        params["conductance"] = Conductance.from_df_state(df, grid_id)
        params["irrigation"] = IrrigationParams.from_df_state(df, grid_id)
        params["anthropogenic_emissions"] = AnthropogenicEmissions.from_df_state(
            df, grid_id
        )
        params["snow"] = SnowParams.from_df_state(df, grid_id)
        params["land_cover"] = LandCover.from_df_state(df, grid_id)
        params["vertical_layers"] = VerticalLayers.from_df_state(df, grid_id)

        params["stebbs"] = StebbsProperties.from_df_state(df, grid_id)
        params["building_archetype"] = ArchetypeProperties.from_df_state(df, grid_id)

        return cls(**params)


class Site(BaseModel):
    name: str = Field(description="Name of the site", default="test site")
    gridiv: int = Field(
        description="Grid ID for identifying this site in multi-site simulations",
        default=1,
    )
    properties: SiteProperties = Field(
        default_factory=SiteProperties,
        description="Physical and morphological properties of the site",
    )
    initial_states: InitialStates = Field(
        default_factory=InitialStates,
        description="Initial conditions for model state variables",
    )

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert site to DataFrame state format"""
        df_state = init_df_state(grid_id)
        df_site_properties = self.properties.to_df_state(grid_id)
        df_initial_states = self.initial_states.to_df_state(grid_id)
        df_state = pd.concat([df_state, df_site_properties, df_initial_states], axis=1)
        return df_state


class SnowAlb(BaseModel):
    snowalb: FlexibleRefValue(float) = Field(
        description="Snow albedo", json_schema_extra={"unit": "dimensionless"},
        default=0.7,
        ge=0,
        le=1,
    )

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert snow albedo to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing snow albedo parameters
        """
        df_state = init_df_state(grid_id)
        df_state[("snowalb", "0")] = self.snowalb.value if isinstance(self.snowalb, RefValue) else self.snowalb
        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "SnowAlb":
        """
        Reconstruct SnowAlb from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing snow albedo parameters.
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            SnowAlb: Instance of SnowAlb.
        """
        snowalb = df.loc[grid_id, ("snowalb", "0")]
        return cls(snowalb=RefValue(snowalb))


class DLSCheck(BaseModel):
    lat: float
    lng: float
    year: int
    startdls: Optional[int] = None
    enddls: Optional[int] = None

    def compute_dst_transitions(self):
        tf = TimezoneFinder()
        tz_name = tf.timezone_at(lat=self.lat, lng=self.lng)

        if not tz_name:
            print(f"[DLS] ❌ Cannot determine timezone for lat={self.lat}, lng={self.lng}")
            return None, None, None

        tz = pytz.timezone(tz_name)

        def find_transition(month: int) -> Optional[int]:
            try:
                prev_dt = tz.localize(datetime(self.year, month, 1, 12), is_dst=None)
                prev_offset = prev_dt.utcoffset()
                for day in range(2, 32):
                    try:
                        curr_dt = tz.localize(datetime(self.year, month, day, 12), is_dst=None)
                        curr_offset = curr_dt.utcoffset()
                        if curr_offset != prev_offset:
                            return curr_dt.timetuple().tm_yday
                        prev_offset = curr_offset
                    except Exception:
                        continue
                return None
            except Exception:
                return None

        if self.lat >= 0:  # Northern Hemisphere
            start = find_transition(3) or find_transition(4)
            end = find_transition(10) or find_transition(11)
        else:  # Southern Hemisphere
            start = find_transition(9) or find_transition(10)
            end = find_transition(3) or find_transition(4)

        return start, end, tz_name


# class DLSValidator(BaseModel):
#     lat: float
#     lng: float
#     startdls: Optional[float] = None
#     enddls: Optional[float] = None
#     year: int 

#     def validate_dls(self):
#         print("CHECKING -- Daylight-Saving dates.")
#         print(f"[DEBUG] Running validate_dls with: self.lat={self.lat}, self.lng={self.lng}, self.startdls={self.startdls}, self.enddls={self.enddls}, self.year={self.year}")

#         tf = TimezoneFinder()
#         tz_name = tf.timezone_at(lat=self.lat, lng=self.lng)
#         if not tz_name:
#             print(f"[DEBUG] Cannot determine timezone for lat={self.lat}, lng={self.lng}; skipping DST check.")
#             return

#         if tz_name.startswith("Etc/"):
#             print(f"[DEBUG] Timezone returned is '{tz_name}', which may not support DST.")
#             return

#         tz = pytz.timezone(tz_name)

#         def find_transition(month: int) -> Optional[int]:
#             try:
#                 prev_dt = tz.localize(datetime(self.year, month, 1, 12, 0), is_dst=None)
#             except Exception:
#                 return None
#             prev_off = prev_dt.utcoffset()

#             for day in range(2, 32):
#                 try:
#                     curr_dt = tz.localize(datetime(self.year, month, day, 12, 0), is_dst=None)
#                 except Exception:
#                     continue
#                 curr_off = curr_dt.utcoffset()
#                 if curr_off != prev_off:
#                     return curr_dt.timetuple().tm_yday
#                 prev_off = curr_off
#             return None

#         if self.lat >= 0:
#             actual_start = find_transition(3) or find_transition(4)
#             actual_end = find_transition(10) or find_transition(11)
#         else:
#             actual_start = find_transition(9) or find_transition(10)
#             actual_end = find_transition(3) or find_transition(4)

#         print(f"[DEBUG] Computed DST start={actual_start}, end={actual_end}")
#         print(f"[DEBUG] Provided startdls={self.startdls}, enddls={self.enddls}")

#         if actual_start and actual_end:
#             # Allow some tolerance for DST dates (within 7 days)
#             start_diff = abs(self.startdls - actual_start) if self.startdls else 0
#             end_diff = abs(self.enddls - actual_end) if self.enddls else 0
#             
#             if start_diff > 7 or end_diff > 7:
#                 raise ValueError(
#                     f"ERROR — DLS mismatch:\n"
#                     f"  computed start={actual_start}, end={actual_end}\n"
#                     f"  in YAML   start={self.startdls}, end={self.enddls}"
#                 )
#             else:
#                 print(f"Daylight-Saving dates OK (tolerance: start_diff={start_diff}, end_diff={end_diff}).")
#         else:
#             print("⚠️ Unable to detect one or both DST transitions; please verify manually.")


class SeasonCheck(BaseModel):
    start_date: str  # YYYY-MM-DD
    end_date: str
    lat: float

    def get_season(self) -> str:
        try:
            start = datetime.strptime(self.start_date, "%Y-%m-%d").timetuple().tm_yday
            end = datetime.strptime(self.end_date, "%Y-%m-%d").timetuple().tm_yday
        except ValueError:
            raise ValueError("start_date and end_date must be in YYYY-MM-DD format")

        abs_lat = abs(self.lat)

        # Near equator: no season
        if abs_lat <= 10:
            return "equatorial"

        # Tropical belt
        if 10 < abs_lat < 23.5:
            return "tropical"

        # Standard seasonal logic
        if self.lat >= 0:  # Northern Hemisphere
            if 150 < start < 250 and 150 < end < 250:
                return "summer"
            elif 60 < start <= 150 and 60 < end <= 150:
                return "spring"
            elif (start <= 60 or start >= 335) and (end <= 60 or end >= 335):
                return "winter"
        else:  # Southern Hemisphere
            if 150 < start < 250 and 150 < end < 250:
                return "winter"
            elif 60 < start <= 150 and 60 < end <= 150:
                return "fall"
            elif 250 <= start < 335 and 250 <= end < 335:
                return "spring"
            elif (start <= 60 or start >= 335) and (end <= 60 or end >= 335):
                return "summer"

        return "unknown"
