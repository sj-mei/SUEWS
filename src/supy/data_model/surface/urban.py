"""Urban surface models for SUEWS."""

from typing import Dict, List, Optional, Union, Literal
import pandas as pd
from pydantic import BaseModel, Field, model_validator, PrivateAttr

from ..base import ValueWithDOI, Reference
from .base import SurfaceType, init_df_state


class WaterDistribution(BaseModel):
    """Water distribution parameters for surfaces."""
    # Optional fields for all possible distributions
    to_paved: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_bldgs: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_dectr: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_evetr: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_grass: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_bsoil: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_water: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_runoff: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)  # For paved/bldgs
    to_soilstore: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)  # For vegetated surfaces
    _surface_type: Optional[SurfaceType] = PrivateAttr(None)

    ref: Optional[Reference] = None

    def __init__(self, surface_type: Optional[SurfaceType] = None, **data):
        # Store surface type as private attribute
        super().__init__(**data)
        self._surface_type = surface_type

        # If surface type is provided, set default values
        if surface_type:
            self._set_defaults(surface_type)
            self.validate_distribution(surface_type)

    def _set_defaults(self, surface_type: SurfaceType):
        # Default distributions based on surface type
        default_distributions = {
            SurfaceType.PAVED: {
                "to_bldgs": ValueWithDOI(0.2),
                "to_evetr": ValueWithDOI(0.1),
                "to_dectr": ValueWithDOI(0.1),
                "to_grass": ValueWithDOI(0.1),
                "to_bsoil": ValueWithDOI(0.1),
                "to_water": ValueWithDOI(0.1),
                "to_runoff": ValueWithDOI(0.3),
            },
            SurfaceType.BLDGS: {
                "to_paved": ValueWithDOI(0.2),
                "to_evetr": ValueWithDOI(0.1),
                "to_dectr": ValueWithDOI(0.1),
                "to_grass": ValueWithDOI(0.1),
                "to_bsoil": ValueWithDOI(0.1),
                "to_water": ValueWithDOI(0.1),
                "to_runoff": ValueWithDOI(0.3),
            },
        }

        if surface_type in default_distributions:
            defaults = default_distributions[surface_type]
            for key, value in defaults.items():
                if getattr(self, key) is None:
                    setattr(self, key, value)

    def validate_distribution(self, surface_type: SurfaceType) -> None:
        """Validate water distribution based on surface type"""
        # Define required fields for each surface type
        required_fields = {
            SurfaceType.PAVED: [
                "to_bldgs",
                "to_dectr",
                "to_evetr",
                "to_grass",
                "to_bsoil",
                "to_water",
                "to_runoff",
            ],
            SurfaceType.BLDGS: [
                "to_paved",
                "to_dectr",
                "to_evetr",
                "to_grass",
                "to_bsoil",
                "to_water",
                "to_runoff",
            ],
        }

        if surface_type == SurfaceType.WATER:
            raise ValueError("Water surface should not have water distribution")

        fields = required_fields[surface_type]
        values = []

        # Check required fields are present and collect values
        for field in fields:
            value = getattr(self, field)
            if value is None:
                raise ValueError(f"Missing required field {field} for {surface_type.value}")
            values.append(value)

        # Validate sum
        total = sum(value.value for value in values)
        if not abs(total - 1.0) < 1e-5:
            raise ValueError(f"Water distribution sum must be 1.0, got {total}")

    def to_df_state(self, grid_id: int, surf_idx: int) -> pd.DataFrame:
        """Convert water distribution parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            surf_idx: Surface index (0=paved, 1=bldgs)

        Returns:
            pd.DataFrame: DataFrame containing water distribution parameters
        """
        df_state = init_df_state(grid_id)

        # Add all non-None distribution parameters using a two-step process
        # 1. Collect all non-None values
        list_waterdist_value = []
        for i, attr in enumerate(
            [
                "to_paved",
                "to_bldgs",
                "to_evetr",
                "to_dectr",
                "to_grass",
                "to_bsoil",
                "to_water",
            ]
        ):
            value = getattr(self, attr)
            if value is None:
                list_waterdist_value.append(0.0)
            else:
                list_waterdist_value.append(value)

        # either to_soilstore or to_runoff must be provided - the other must be 0
        to_soilstore_or_runoff = (
            self.to_runoff if self.to_soilstore is None else self.to_soilstore
        )
        list_waterdist_value.append(to_soilstore_or_runoff)

        # 2. Create param_tuples and values - only add non-None values following the order of the list
        for i, value in enumerate(list_waterdist_value):
            if value is not None:
                df_state.loc[grid_id, ("waterdist", f"({i}, {surf_idx})")] = (
                    value.value if isinstance(value, ValueWithDOI) else value
                )

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int, surf_idx: int) -> "WaterDistribution":
        """Reconstruct WaterDistribution from DataFrame state format.

        Args:
            df: DataFrame containing water distribution parameters
            grid_id: Grid ID for the DataFrame index
            surf_idx: Surface index for identifying columns

        Returns:
            WaterDistribution: Instance of WaterDistribution
        """
        dict_surface_type = {
            0: SurfaceType.PAVED,
            1: SurfaceType.BLDGS,
        }
        surface_type = dict_surface_type[surf_idx]
        # initialize an instance of this class
        instance = cls(surface_type=surface_type)

        # Define the parameter names and their indices
        param_map = {
            "to_paved": 0,
            "to_bldgs": 1,
            "to_evetr": 2,
            "to_dectr": 3,
            "to_grass": 4,
            "to_bsoil": 5,
            "to_water": 6,
        }

        # Extract the values from the DataFrame
        params = {
            param: df.loc[grid_id, ("waterdist", f"({idx}, {surf_idx})")]
            for param, idx in param_map.items()
        }
        for param, value in params.items():
            value = ValueWithDOI(value)
            if getattr(instance, param) is not None:
                setattr(instance, param, value)

        # set the last to_soilstore or to_runoff
        waterdist_last = df.loc[grid_id, ("waterdist", f"(7, {surf_idx})")]
        waterdist_last = ValueWithDOI(waterdist_last)
        if getattr(instance, "to_soilstore") is None:
            setattr(instance, "to_runoff", waterdist_last)
        else:
            setattr(instance, "to_soilstore", waterdist_last)

        return instance


class StorageDrainParams(BaseModel):
    """Storage and drainage parameters for surfaces."""
    store_min: ValueWithDOI[float] = Field(ge=0, default=ValueWithDOI(0.0))
    store_max: ValueWithDOI[float] = Field(ge=0, default=ValueWithDOI(10.0))
    store_cap: ValueWithDOI[float] = Field(ge=0, default=ValueWithDOI(10.0))
    drain_eq: ValueWithDOI[int] = Field(default=ValueWithDOI(0))
    drain_coef_1: ValueWithDOI[float] = Field(default=ValueWithDOI(0.013))
    drain_coef_2: ValueWithDOI[float] = Field(default=ValueWithDOI(1.71))

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int, surf_idx: int) -> pd.DataFrame:
        """Convert storage and drain parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            surf_idx: Surface index (0=paved, 1=bldgs)

        Returns:
            pd.DataFrame: DataFrame containing storage and drain parameters
        """
        df_state = init_df_state(grid_id)

        # Set parameters in order
        for i, value in enumerate(
            [
                self.store_min,
                self.drain_eq,
                self.drain_coef_1,
                self.drain_coef_2,
                self.store_max,
                self.store_cap,
            ]
        ):
            df_state.loc[grid_id, ("storedrainprm", f"({i}, {surf_idx})")] = value.value

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int, surf_idx: int) -> "StorageDrainParams":
        """Reconstruct StorageDrainParams from DataFrame state format.

        Args:
            df: DataFrame containing storage and drain parameters
            grid_id: Grid ID for the DataFrame index
            surf_idx: Surface index for identifying columns

        Returns:
            StorageDrainParams: Instance of StorageDrainParams
        """
        # Map each parameter to its corresponding index
        param_map = {
            "store_min": 0,
            "drain_eq": 1,
            "drain_coef_1": 2,
            "drain_coef_2": 3,
            "store_max": 4,
            "store_cap": 5,
        }

        # Extract values from DataFrame
        params = {
            param: df.loc[grid_id, ("storedrainprm", f"({idx}, {surf_idx})")]
            for param, idx in param_map.items()
        }

        # Convert to ValueWithDOI
        params = {key: ValueWithDOI(value) for key, value in params.items()}

        return cls(**params)


class NonVegetatedSurfaceProperties(BaseModel):
    """Base properties for non-vegetated surfaces."""
    sfr: ValueWithDOI[float] = Field(
        ge=0, le=1,
        description="Surface fraction",
        default=ValueWithDOI(1.0 / 7),
    )
    emis: ValueWithDOI[float] = Field(
        ge=0, le=1,
        description="Surface emissivity",
        default=ValueWithDOI(0.95),
    )
    chanohm: Optional[ValueWithDOI[float]] = Field(default=ValueWithDOI(0.0))
    cpanohm: Optional[ValueWithDOI[float]] = Field(default=ValueWithDOI(1200.0))
    kkanohm: Optional[ValueWithDOI[float]] = Field(default=ValueWithDOI(0.4))
    ohm_threshsw: Optional[ValueWithDOI[float]] = Field(default=ValueWithDOI(0.0))
    ohm_threshwd: Optional[ValueWithDOI[float]] = Field(default=ValueWithDOI(0.0))
    soildepth: ValueWithDOI[float] = Field(default=ValueWithDOI(0.15))
    soilstorecap: ValueWithDOI[float] = Field(default=ValueWithDOI(150.0))
    statelimit: ValueWithDOI[float] = Field(default=ValueWithDOI(10.0))
    wetthresh: ValueWithDOI[float] = Field(default=ValueWithDOI(0.5))
    sathydraulicconduct: ValueWithDOI[float] = Field(default=ValueWithDOI(0.0001))
    waterdist: Optional[WaterDistribution] = Field(
        default=None,
        description="Water distribution parameters",
    )
    storedrainprm: StorageDrainParams = Field(
        default_factory=StorageDrainParams,
        description="Storage and drain parameters",
    )
    snowpacklimit: Optional[ValueWithDOI[float]] = Field(default=ValueWithDOI(10.0))
    irrfrac: Optional[ValueWithDOI[float]] = Field(default=ValueWithDOI(0.0))
    _surface_type: Optional[SurfaceType] = PrivateAttr(default=None)

    ref: Optional[Reference] = None

    def set_surface_type(self, surface_type: SurfaceType):
        self._surface_type = surface_type
        if self._surface_type == SurfaceType.WATER:
            if self.waterdist is not None:
                raise ValueError("Water surface should not have water distribution")
        else:
            if self.waterdist is None:
                raise ValueError(f"Water distribution required for {self._surface_type.value}")
            self.waterdist.validate_distribution(self._surface_type)

    def get_surface_type(self) -> SurfaceType:
        return self._surface_type

    def get_surface_name(self) -> str:
        return self._surface_type.value

    def get_surface_index(self) -> int:
        dict_surface_type = {
            SurfaceType.PAVED: 0,
            SurfaceType.BLDGS: 1,
        }
        return dict_surface_type[self._surface_type]

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert surface properties to DataFrame state format."""
        df_state = init_df_state(grid_id)

        # Get surface index
        surf_idx = self.get_surface_index()

        # Get surface name
        surf_name = self.get_surface_name()

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = f"({surf_idx},)"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = None
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # Add all non-inherited properties
        for attr in [
            "sfr",
            "emis",
            "chanohm",
            "cpanohm",
            "kkanohm",
            "ohm_threshsw",
            "ohm_threshwd",
            "soildepth",
            "soilstorecap",
            "statelimit",
            "wetthresh",
            "sathydraulicconduct",
            "snowpacklimit",
        ]:
            value = getattr(self, attr)
            if value is not None:
                set_df_value(attr, value.value)

        # Add water distribution if present
        if self.waterdist is not None:
            df_waterdist = self.waterdist.to_df_state(grid_id, surf_idx)
            df_state = pd.concat([df_state, df_waterdist], axis=1)

        # Add storage drain parameters
        df_storedrainprm = self.storedrainprm.to_df_state(grid_id, surf_idx)
        df_state = pd.concat([df_state, df_storedrainprm], axis=1)

        # Add irrfrac
        df_state.loc[grid_id, (f"irrfrac{surf_name}", "0")] = self.irrfrac.value

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int, surf_idx: int) -> "NonVegetatedSurfaceProperties":
        """Reconstruct surface properties from DataFrame state format."""
        instance = cls()

        # Set surface type based on index
        surface_type = {0: SurfaceType.PAVED, 1: SurfaceType.BLDGS}[surf_idx]
        instance.set_surface_type(surface_type)

        # Helper function to get values from DataFrame
        def get_df_value(col_name: str) -> float:
            return df.loc[grid_id, (col_name, f"({surf_idx},)")]

        # Set all non-inherited properties
        for attr in [
            "sfr",
            "emis",
            "chanohm",
            "cpanohm",
            "kkanohm",
            "ohm_threshsw",
            "ohm_threshwd",
            "soildepth",
            "soilstorecap",
            "statelimit",
            "wetthresh",
            "sathydraulicconduct",
            "snowpacklimit",
        ]:
            try:
                value = get_df_value(attr)
                setattr(instance, attr, ValueWithDOI(value))
            except KeyError:
                continue

        # Set water distribution
        instance.waterdist = WaterDistribution.from_df_state(df, grid_id, surf_idx)

        # Set storage drain parameters
        instance.storedrainprm = StorageDrainParams.from_df_state(df, grid_id, surf_idx)

        # Set irrfrac
        surf_name = instance.get_surface_name()
        instance.irrfrac = ValueWithDOI(df.loc[grid_id, (f"irrfrac{surf_name}", "0")])

        return instance


class PavedProperties(NonVegetatedSurfaceProperties):
    """Properties for paved surfaces."""
    _surface_type: Literal[SurfaceType.PAVED] = SurfaceType.PAVED
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.PAVED),
        description="Water distribution for paved surfaces",
    )

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert paved surface properties to DataFrame state format."""
        return super().to_df_state(grid_id)

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "PavedProperties":
        """Reconstruct paved surface properties from DataFrame state format."""
        surf_idx = 0
        return super().from_df_state(df, grid_id, surf_idx)


class BldgsProperties(NonVegetatedSurfaceProperties):
    """Properties for building surfaces."""
    _surface_type: Literal[SurfaceType.BLDGS] = SurfaceType.BLDGS
    faibldg: ValueWithDOI[float] = Field(
        ge=0,
        default=ValueWithDOI(0.3),
        description="Frontal area index of buildings",
    )
    bldgh: ValueWithDOI[float] = Field(
        ge=3,
        default=ValueWithDOI(10.0),
        description="Building height",
    )
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.BLDGS),
        description="Water distribution for buildings",
    )

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert building properties to DataFrame state format."""
        df_state = super().to_df_state(grid_id)

        # Add building-specific properties
        df_state.loc[grid_id, ("faibldg", "0")] = self.faibldg.value
        df_state.loc[grid_id, ("bldgh", "0")] = self.bldgh.value

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "BldgsProperties":
        """Reconstruct building properties from DataFrame state format."""
        surf_idx = 1
        instance = super().from_df_state(df, grid_id, surf_idx)

        # Set building-specific properties
        instance.faibldg = ValueWithDOI(df.loc[grid_id, ("faibldg", "0")])
        instance.bldgh = ValueWithDOI(df.loc[grid_id, ("bldgh", "0")])

        return instance
