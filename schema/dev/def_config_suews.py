from typing import Dict, List, Optional, Union, Literal, Tuple
from pydantic import (
    BaseModel,
    Field,
    model_validator,
    field_validator,
    PrivateAttr,
)
import numpy as np
from enum import Enum
import pandas as pd
import scipy as sp
import yaml
import pdb
import math


def init_df_state(grid_id: int) -> pd.DataFrame:
    idx = pd.Index([grid_id], name="grid")
    col = pd.MultiIndex.from_tuples([("grid_iv", 0)], names=["var", "ind_dim"])
    df_state = pd.DataFrame(index=idx, columns=col)
    df_state.loc[grid_id, ("grid_iv", 0)] = grid_id
    return df_state


class SurfaceType(str, Enum):
    PAVED = "paved"
    BLDGS = "bldgs"
    EVETR = "evetr"
    DECTR = "dectr"
    GRASS = "grass"
    BSOIL = "bsoil"
    WATER = "water"


class SnowAlb(BaseModel):
    snowalb: float = Field(ge=0, le=1, description="Snow albedo", default=0.5)

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert snow albedo to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing snow albedo parameters
        """
        df_state = init_df_state(grid_id)
        df_state[("snowalb", "0")] = self.snowalb
        return df_state


class WaterUse(BaseModel):
    wu_total: float = Field(ge=0, description="Total water use", default=0.0)
    wu_auto: float = Field(ge=0, description="Automatic water use", default=0.0)
    wu_manual: float = Field(ge=0, description="Manual water use", default=0.0)

    def to_df_state(self, veg_idx: int, grid_id: int) -> pd.DataFrame:
        """Convert water use to DataFrame state format."""
        df_state = init_df_state(grid_id)
        df_state.loc[grid_id, ("wuday_id", f"{veg_idx * 3 + 0}")] = self.wu_total
        df_state.loc[grid_id, ("wuday_id", f"{veg_idx * 3 + 1}")] = self.wu_auto
        df_state.loc[grid_id, ("wuday_id", f"{veg_idx * 3 + 2}")] = self.wu_manual
        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, veg_idx: int, grid_id: int) -> "WaterUse":
        """
        Reconstruct WaterUse from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing water use parameters.
            veg_idx (int): Vegetation index for identifying columns.
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            WaterUse: Instance of WaterUse.
        """
        wu_total = df.loc[grid_id, ("wuday_id", f"{veg_idx * 3 + 0}")]
        wu_auto = df.loc[grid_id, ("wuday_id", f"{veg_idx * 3 + 1}")]
        wu_manual = df.loc[grid_id, ("wuday_id", f"{veg_idx * 3 + 2}")]

        return cls(wu_total=wu_total, wu_auto=wu_auto, wu_manual=wu_manual)


class SurfaceInitialState(BaseModel):
    """Base initial state parameters for all surface types"""

    state: float = Field(ge=0, description="Initial state of the surface", default=0.0)
    soilstore: float = Field(ge=0, description="Initial soil store", default=0.0)
    snowfrac: Optional[float] = Field(
        ge=0, le=1, description="Snow fraction", default=0.0
    )
    snowpack: Optional[float] = Field(ge=0, description="Snow pack", default=0.0)
    icefrac: Optional[float] = Field(
        ge=0, le=1, description="Ice fraction", default=0.0
    )
    snowwater: Optional[float] = Field(ge=0, description="Snow water", default=0.0)
    snowdens: Optional[float] = Field(ge=0, description="Snow density", default=0.0)
    temperature: List[float] = Field(
        min_items=5,
        max_items=5,
        description="Initial temperature for each thermal layer",
        default=[15.0, 15.0, 15.0, 15.0, 15.0],
    )
    tsfc: Optional[float] = Field(
        description="Initial exterior surface temperature", default=15.0
    )
    tin: Optional[float] = Field(
        description="Initial interior surface temperature", default=20.0
    )
    _surface_type: Optional[SurfaceType] = PrivateAttr(default=None)

    def set_surface_type(self, surface_type: SurfaceType):
        """Set surface type"""
        self._surface_type = surface_type

    def get_surface_index(self) -> int:
        """Get surface index"""
        return {
            SurfaceType.PAVED: 0,
            SurfaceType.BLDGS: 1,
            SurfaceType.EVETR: 2,
            SurfaceType.DECTR: 3,
            SurfaceType.GRASS: 4,
            SurfaceType.BSOIL: 5,
            SurfaceType.WATER: 6,
        }[self._surface_type]

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert base surface initial state to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing initial state parameters
        """
        df_state = init_df_state(grid_id)

        # Get surface index
        surf_idx = self.get_surface_index()

        # Set basic state parameters
        df_state[("state", f"({surf_idx},)")] = self.state
        df_state[("soilstore_id", f"({surf_idx},)")] = self.soilstore

        # Set snow/ice parameters if present
        if self.snowfrac is not None:
            df_state[("snowfrac", f"({surf_idx},)")] = self.snowfrac
        if self.snowpack is not None:
            df_state[("snowpack", f"({surf_idx},)")] = self.snowpack
        if self.icefrac is not None:
            df_state[("icefrac", f"({surf_idx},)")] = self.icefrac
        if self.snowwater is not None:
            df_state[("snowwater", f"({surf_idx},)")] = self.snowwater
        if self.snowdens is not None:
            df_state[("snowdens", f"({surf_idx},)")] = self.snowdens

        # Set temperature parameters
        for i, temp in enumerate(self.temperature):
            df_state[("t", f"({surf_idx},{i})")] = temp

        if self.tsfc is not None:
            df_state[("tsfc", f"({surf_idx},)")] = self.tsfc
        if self.tin is not None:
            df_state[("tin", f"({surf_idx},)")] = self.tin

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int
    ) -> "SurfaceInitialState":
        """
        Reconstruct SurfaceInitialState from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing surface initial state parameters.
            grid_id (int): Grid ID for the DataFrame index.
            surf_idx (int): Surface index for identifying columns.

        Returns:
            SurfaceInitialState: Instance of SurfaceInitialState.
        """
        state = df.loc[grid_id, ("state", f"({surf_idx},)")]
        soilstore = df.loc[grid_id, ("soilstore_id", f"({surf_idx},)")]
        snowfrac = df.loc[grid_id, ("snowfrac", f"({surf_idx},)")]
        snowpack = df.loc[grid_id, ("snowpack", f"({surf_idx},)")]
        icefrac = df.loc[grid_id, ("icefrac", f"({surf_idx},)")]
        snowwater = df.loc[grid_id, ("snowwater", f"({surf_idx},)")]
        snowdens = df.loc[grid_id, ("snowdens", f"({surf_idx},)")]

        temperature = [df.loc[grid_id, ("t", f"({surf_idx},{i})")] for i in range(5)]
        tsfc = df.loc[grid_id, ("tsfc", f"({surf_idx},)")]
        tin = df.loc[grid_id, ("tin", f"({surf_idx},)")]

        return cls(
            state=state,
            soilstore=soilstore,
            snowfrac=snowfrac,
            snowpack=snowpack,
            icefrac=icefrac,
            snowwater=snowwater,
            snowdens=snowdens,
            temperature=temperature,
            tsfc=tsfc,
            tin=tin,
        )


class VegetatedSurfaceInitialState(SurfaceInitialState):
    """Base initial state parameters for vegetated surfaces"""

    alb_id: float = Field(
        description="Initial albedo for vegetated surfaces", default=0.1
    )
    lai_id: float = Field(description="Initial leaf area index", default=1.0)
    gdd_id: float = Field(description="Growing degree days ID", default=0)
    sdd_id: float = Field(description="Senescence degree days ID", default=0)
    wu: WaterUse = Field(default_factory=WaterUse)

    @model_validator(mode="after")
    def validate_surface_state(self) -> "VegetatedSurfaceInitialState":
        """Validate state based on surface type"""
        # Skip validation if surface type not yet set
        if not hasattr(self, "_surface_type") or self._surface_type is None:
            return self

        if self._surface_type not in [
            SurfaceType.DECTR,
            SurfaceType.EVETR,
            SurfaceType.GRASS,
        ]:
            raise ValueError(
                f"Invalid surface type {self._surface_type} for vegetated surface"
            )
        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert vegetated surface initial state to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing initial state parameters
        """
        # Get base surface state parameters
        df_state = super().to_df_state(grid_id)

        # Get surface index
        surf_idx = self.get_surface_index()

        # Add vegetated surface specific parameters
        df_state[("alb", f"({surf_idx},)")] = self.alb_id
        df_state[("lai", f"({surf_idx},)")] = self.lai_id
        df_state[("gdd", f"({surf_idx},)")] = self.gdd_id
        df_state[("sdd", f"({surf_idx},)")] = self.sdd_id

        # Add water use parameters
        veg_idx = surf_idx - 2
        df_wu = self.wu.to_df_state(veg_idx, grid_id)
        df_state = pd.concat([df_state, df_wu], axis=1)

        # Drop any duplicate columns
        df_state = df_state.loc[:, ~df_state.columns.duplicated()]

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int
    ) -> "VegetatedSurfaceInitialState":
        """
        Reconstruct VegetatedSurfaceInitialState from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing vegetated surface state parameters.
            grid_id (int): Grid ID for the DataFrame index.
            surf_idx (int): Surface index for identifying columns.

        Returns:
            VegetatedSurfaceInitialState: Instance of VegetatedSurfaceInitialState.
        """
        # Base class reconstruction
        base_instance = SurfaceInitialState.from_df_state(df, grid_id, surf_idx)

        # Vegetated surface-specific parameters
        alb_id = df.loc[grid_id, ("alb", f"({surf_idx},)")]
        lai_id = df.loc[grid_id, ("lai", f"({surf_idx},)")]
        gdd_id = df.loc[grid_id, ("gdd", f"({surf_idx},)")]
        sdd_id = df.loc[grid_id, ("sdd", f"({surf_idx},)")]

        # Reconstruct WaterUse instance
        veg_idx = surf_idx - 2
        wu = WaterUse.from_df_state(df, veg_idx, grid_id)

        return cls(
            **base_instance.dict(),
            alb_id=alb_id,
            lai_id=lai_id,
            gdd_id=gdd_id,
            sdd_id=sdd_id,
            wu=wu,
        )


class DeciduousTreeSurfaceInitialState(VegetatedSurfaceInitialState):
    """Initial state parameters for deciduous trees"""

    porosity_id: float = Field(description="Initial porosity for deciduous trees")
    decidcap_id: float = Field(
        description="Initial deciduous capacity for deciduous trees"
    )

    @model_validator(mode="after")
    def validate_surface_state(self) -> "DeciduousTreeSurfaceInitialState":
        """Validate state based on surface type"""
        # Skip validation if surface type not yet set
        if not hasattr(self, "_surface_type") or self._surface_type is None:
            return self

        if self._surface_type != SurfaceType.DECTR:
            raise ValueError("This state is only valid for deciduous trees")
        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert deciduous tree initial state to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing initial state parameters
        """
        # Get base vegetated surface state parameters
        df_state = super().to_df_state(grid_id)

        # Get surface index
        surf_idx = self.get_surface_index()

        # Add deciduous tree specific parameters
        df_state[("porosity", f"({surf_idx},)")] = self.porosity_id
        df_state[("decidcap", f"({surf_idx},)")] = self.decidcap_id

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int
    ) -> "DeciduousTreeSurfaceInitialState":
        """
        Reconstruct DeciduousTreeSurfaceInitialState from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing deciduous tree state parameters.
            grid_id (int): Grid ID for the DataFrame index.
            surf_idx (int): Surface index for identifying columns.

        Returns:
            DeciduousTreeSurfaceInitialState: Instance of DeciduousTreeSurfaceInitialState.
        """
        # Base class reconstruction
        base_instance = VegetatedSurfaceInitialState.from_df_state(
            df, grid_id, surf_idx
        )

        # Deciduous tree-specific parameters
        porosity_id = df.loc[grid_id, ("porosity", f"({surf_idx},)")]
        decidcap_id = df.loc[grid_id, ("decidcap", f"({surf_idx},)")]

        return cls(
            **base_instance.dict(),
            porosity_id=porosity_id,
            decidcap_id=decidcap_id,
        )


class InitialStates(BaseModel):
    """Initial conditions for the SUEWS model"""

    snowalb: float = Field(ge=0, le=1, description="Initial snow albedo", default=0.5)
    paved: SurfaceInitialState = Field(default_factory=SurfaceInitialState)
    bldgs: SurfaceInitialState = Field(default_factory=SurfaceInitialState)
    evetr: VegetatedSurfaceInitialState = Field(
        default_factory=VegetatedSurfaceInitialState
    )
    dectr: DeciduousTreeSurfaceInitialState = Field(
        default_factory=DeciduousTreeSurfaceInitialState
    )
    grass: VegetatedSurfaceInitialState = Field(
        default_factory=VegetatedSurfaceInitialState
    )
    bsoil: SurfaceInitialState = Field(default_factory=SurfaceInitialState)
    water: SurfaceInitialState = Field(default_factory=SurfaceInitialState)
    roofs: Optional[List[SurfaceInitialState]] = Field(
        default=[
            SurfaceInitialState(),
            SurfaceInitialState(),
        ],
        description="Initial states for roof layers",
    )
    walls: Optional[List[SurfaceInitialState]] = Field(
        default=[
            SurfaceInitialState(),
            SurfaceInitialState(),
        ],
        description="Initial states for wall layers",
    )

    @model_validator(mode="before")
    @classmethod
    def set_surface_types(cls, data: Dict) -> Dict:
        """Set surface types for all surfaces before validation"""
        # Create instances if they don't exist
        for surface_type in ["paved", "bldgs", "bsoil", "water"]:
            if surface_type not in data:
                data[surface_type] = SurfaceInitialState()
            if isinstance(data[surface_type], dict):
                data[surface_type] = SurfaceInitialState(**data[surface_type])
            data[surface_type].set_surface_type(SurfaceType(surface_type))

        # Handle vegetated surfaces
        if "evetr" in data:
            if isinstance(data["evetr"], dict):
                data["evetr"] = VegetatedSurfaceInitialState(**data["evetr"])
            data["evetr"].set_surface_type(SurfaceType.EVETR)

        if "dectr" in data:
            if isinstance(data["dectr"], dict):
                data["dectr"] = DeciduousTreeSurfaceInitialState(**data["dectr"])
            data["dectr"].set_surface_type(SurfaceType.DECTR)

        if "grass" in data:
            if isinstance(data["grass"], dict):
                data["grass"] = VegetatedSurfaceInitialState(**data["grass"])
            data["grass"].set_surface_type(SurfaceType.GRASS)

        return data

    # def __init__(self, **data):
    #     super().__init__(**data)
    #     # Set surface types for non-vegetated surfaces
    #     self.paved.set_surface_type(SurfaceType.PAVED)
    #     self.bldgs.set_surface_type(SurfaceType.BLDGS)
    #     self.bsoil.set_surface_type(SurfaceType.BSOIL)
    #     self.water.set_surface_type(SurfaceType.WATER)

    #     # Set surface types for vegetated surfaces
    #     # These need to be set before validation
    #     if isinstance(self.evetr, VegetatedSurfaceInitialState):
    #         self.evetr.set_surface_type(SurfaceType.EVETR)
    #     if isinstance(self.dectr, DeciduousTreeSurfaceInitialState):
    #         self.dectr.set_surface_type(SurfaceType.DECTR)
    #     if isinstance(self.grass, VegetatedSurfaceInitialState):
    #         self.grass.set_surface_type(SurfaceType.GRASS)

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert initial states to DataFrame state format."""
        df_state = init_df_state(grid_id)

        # Add snowalb
        df_state[("snowalb", "0")] = self.snowalb

        # Add surface states
        surfaces = {
            "paved": self.paved,
            "bldgs": self.bldgs,
            "evetr": self.evetr,
            "dectr": self.dectr,
            "grass": self.grass,
            "bsoil": self.bsoil,
            "water": self.water,
        }
        # Add surface states
        for surface in surfaces.values():
            df_surface = surface.to_df_state(grid_id)
            df_state = pd.concat([df_state, df_surface], axis=1)

        # Add roof and wall states
        for surface_list, surface_type in [(self.roofs, "roof"), (self.walls, "wall")]:
            if surface_list is not None:  # Check for None explicitly
                for i, surface in enumerate(surface_list):
                    if surface is not None:  # Check each surface is not None
                        df_surface = surface.to_df_state(grid_id)
                        # Prefix column names with surface type
                        df_surface.columns = pd.MultiIndex.from_tuples(
                            [
                                (f"{surface_type}_{col[0]}", f"({i},)")
                                for col in df_surface.columns
                            ],
                            names=["var", "ind_dim"],
                        )
                        df_state = pd.concat([df_state, df_surface], axis=1)

        # Drop duplicate columns while preserving first occurrence
        df_state = df_state.loc[:, ~df_state.columns.duplicated(keep="first")]

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "InitialStates":
        """
        Reconstruct InitialStates from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing initial states.
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            InitialStates: Instance of InitialStates.
        """
        snowalb = df.loc[grid_id, ("snowalb", "0")]

        # Reconstruct surface states
        surface_types = {
            "paved": SurfaceInitialState,
            "bldgs": SurfaceInitialState,
            "evetr": VegetatedSurfaceInitialState,
            "dectr": DeciduousTreeSurfaceInitialState,
            "grass": VegetatedSurfaceInitialState,
            "bsoil": SurfaceInitialState,
            "water": SurfaceInitialState,
        }
        surfaces = {
            name: surface_class.from_df_state(df, grid_id, idx)
            for idx, (name, surface_class) in enumerate(surface_types.items())
        }

        # Reconstruct roof and wall states
        def reconstruct_layers(
            layer_name: str, surface_class: Type[SurfaceInitialState]
        ):
            layers = []
            for i in range(2):  # Assuming two layers for simplicity
                try:
                    layer = surface_class.from_df_state(df, grid_id, i)
                    layers.append(layer)
                except KeyError:
                    break
            return layers

        roofs = reconstruct_layers("roof", SurfaceInitialState)
        walls = reconstruct_layers("wall", SurfaceInitialState)

        return cls(
            snowalb=snowalb,
            roofs=roofs,
            walls=walls,
            **surfaces,
        )


class ThermalLayers(BaseModel):
    dz: List[float] = Field([0.1, 0.2, 0.3, 0.4, 0.5], min_items=5, max_items=5)
    k: List[float] = Field([1.0, 1.0, 1.0, 1.0, 1.0], min_items=5, max_items=5)
    cp: List[float] = Field([1000, 1000, 1000, 1000, 1000], min_items=5, max_items=5)

    def to_df_state(self, grid_id: int, surf_idx: int) -> pd.DataFrame:
        """Convert thermal layer parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            surf_idx: Surface index for identifying columns

        Returns:
            pd.DataFrame: DataFrame containing thermal layer parameters
        """
        df_state = init_df_state(grid_id)

        # Add thermal layer parameters
        for i in range(5):
            df_state[("dz", f"({surf_idx},{i})")] = self.dz[i]
            df_state[("k", f"({surf_idx},{i})")] = self.k[i]
            df_state[("cp", f"({surf_idx},{i})")] = self.cp[i]

        return df_state


class VegetationParams(BaseModel):
    porosity_id: int
    gdd_id: int = Field(description="Growing degree days ID")
    sdd_id: int = Field(description="Senescence degree days ID")
    lai: Dict[str, Union[float, List[float]]] = Field(
        description="Leaf area index parameters"
    )
    ie_a: float = Field(description="Irrigation efficiency coefficient a")
    ie_m: float = Field(description="Irrigation efficiency coefficient m")


class WaterDistribution(BaseModel):
    # Optional fields for all possible distributions
    to_paved: Optional[float] = Field(None, ge=0, le=1)
    to_bldgs: Optional[float] = Field(None, ge=0, le=1)
    to_dectr: Optional[float] = Field(None, ge=0, le=1)
    to_evetr: Optional[float] = Field(None, ge=0, le=1)
    to_grass: Optional[float] = Field(None, ge=0, le=1)
    to_bsoil: Optional[float] = Field(None, ge=0, le=1)
    to_water: Optional[float] = Field(None, ge=0, le=1)
    to_runoff: Optional[float] = Field(None, ge=0, le=1)  # For paved/bldgs
    to_soilstore: Optional[float] = Field(None, ge=0, le=1)  # For vegetated surfaces

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
            SurfaceType.DECTR: [
                "to_paved",
                "to_bldgs",
                "to_evetr",
                "to_grass",
                "to_bsoil",
                "to_water",
                "to_soilstore",
            ],
            SurfaceType.EVETR: [
                "to_paved",
                "to_bldgs",
                "to_dectr",
                "to_grass",
                "to_bsoil",
                "to_water",
                "to_soilstore",
            ],
            SurfaceType.GRASS: [
                "to_paved",
                "to_bldgs",
                "to_dectr",
                "to_evetr",
                "to_bsoil",
                "to_water",
                "to_soilstore",
            ],
            SurfaceType.BSOIL: [
                "to_paved",
                "to_bldgs",
                "to_dectr",
                "to_evetr",
                "to_grass",
                "to_water",
                "to_soilstore",
            ],
            SurfaceType.WATER: None,  # Water surface doesn't have water distribution
        }

        if surface_type == SurfaceType.WATER:
            raise ValueError("Water surface should not have water distribution")

        fields = required_fields[surface_type]
        values = []

        # Check required fields are present and collect values
        for field in fields:
            value = getattr(self, field)
            if value is None:
                raise ValueError(
                    f"Missing required field {field} for {surface_type.value}"
                )
            values.append(value)

        # Validate sum
        total = sum(values)
        if not np.isclose(total, 1.0, rtol=1e-5):
            raise ValueError(f"Water distribution sum must be 1.0, got {total}")

    def to_df_state(self, grid_id: int, surf_idx: int) -> pd.DataFrame:
        """Convert water distribution parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            surf_idx: Surface index (0=paved, 1=bldgs, 2=dectr, 3=evetr, 4=grass, 5=bsoil, 6=water)

        Returns:
            pd.DataFrame: DataFrame containing water distribution parameters with MultiIndex columns
        """
        # Create tuples for MultiIndex columns
        param_tuples = []
        values = []

        # Add all non-None distribution parameters
        for i, attr in enumerate(
            [
                "to_paved",
                "to_bldgs",
                "to_dectr",
                "to_evetr",
                "to_grass",
                "to_bsoil",
                "to_water",
                "to_runoff",
                "to_soilstore",
            ]
        ):
            value = getattr(self, attr)
            if value is not None:
                param_tuples.append(("waterdist", (i, surf_idx)))
                values.append(value)

        # Create MultiIndex columns
        columns = pd.MultiIndex.from_tuples(param_tuples, names=["var", "ind_dim"])

        # Create DataFrame with single row
        df = pd.DataFrame(
            index=pd.Index([grid_id], name="grid"),
            columns=columns,
            data=[values],
            dtype=float,
        )

        return df


class StorageDrainParams(BaseModel):
    store_min: float = Field(ge=0, default=0.0)
    store_max: float = Field(ge=0, default=10.0)
    store_cap: float = Field(ge=0, default=10.0)
    drain_eq: int = Field(default=0)
    drain_coef_1: float = Field(default=0.013)
    drain_coef_2: float = Field(default=1.71)

    def to_df_state(self, grid_id: int, surf_idx: int) -> pd.DataFrame:
        """Convert storage and drain parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            surf_idx: Surface index (0=paved, 1=bldgs, 2=dectr, 3=evetr, 4=grass, 5=bsoil, 6=water)

        Returns:
            pd.DataFrame: DataFrame containing storage and drain parameters with MultiIndex columns
        """
        # Create tuples for MultiIndex columns
        param_tuples = [
            ("storedrainprm", f"({i},{surf_idx})")
            for i, _ in enumerate(
                [
                    "store_min",
                    "store_max",
                    "store_cap",
                    "drain_eq",
                    "drain_coef_1",
                    "drain_coef_2",
                ]
            )
        ]

        # Create MultiIndex columns
        columns = pd.MultiIndex.from_tuples(param_tuples, names=["var", "ind_dim"])

        # Create DataFrame with single row
        df = pd.DataFrame(
            index=pd.Index([grid_id], name="grid"), columns=columns, dtype=float
        )

        # Fill values
        for i, var in enumerate(
            [
                "store_min",
                "store_max",
                "store_cap",
                "drain_eq",
                "drain_coef_1",
                "drain_coef_2",
            ]
        ):
            df.loc[grid_id, ("storedrainprm", f"({i},{surf_idx})")] = getattr(self, var)

        return df


class OHM_Coefficient_season_wetness(BaseModel):
    summer_dry: float
    summer_wet: float
    winter_dry: float
    winter_wet: float

    def to_df_state(self, grid_id: int, surf_idx: int, idx_a: int) -> pd.DataFrame:
        """Convert OHM coefficients to DataFrame state format.

        Args:
            grid_id (int): Grid ID
            surf_idx (int): Surface index
            idx_a (int): Index for coefficient (0=a1, 1=a2, 2=a3)

        Returns:
            pd.DataFrame: DataFrame containing OHM coefficients with MultiIndex columns
        """
        df_state = init_df_state(grid_id)

        # Map season/wetness combinations to indices
        season_wetdry_map = {
            "summer_dry": 0,
            "summer_wet": 1,
            "winter_dry": 2,
            "winter_wet": 3,
        }

        # Set values for each season/wetness combination
        for season_wetdry, idx in season_wetdry_map.items():
            str_idx = f"({surf_idx},{idx},{idx_a})"
            df_state.loc[grid_id, ("ohm_coef", str_idx)] = getattr(self, season_wetdry)

        return df_state


class OHMCoefficients(BaseModel):
    a1: OHM_Coefficient_season_wetness
    a2: OHM_Coefficient_season_wetness
    a3: OHM_Coefficient_season_wetness

    def to_df_state(self, grid_id: int, surf_idx: int) -> pd.DataFrame:
        """Convert OHM coefficients to DataFrame state format.

        Args:
            grid_id (int): Grid ID
            surf_idx (int): Surface index

        Returns:
            pd.DataFrame: DataFrame containing OHM coefficients with MultiIndex columns
        """
        df_state = init_df_state(grid_id)

        # Convert each coefficient (a1, a2, a3)
        for idx_a, coef in enumerate([self.a1, self.a2, self.a3]):
            df_coef = coef.to_df_state(grid_id, surf_idx, idx_a)
            df_state = pd.concat([df_state, df_coef], axis=1)

        # drop duplicate columns
        df_state = df_state.loc[:, ~df_state.columns.duplicated()]

        return df_state


class SurfaceProperties(BaseModel):
    """Base properties for all surface types"""

    sfr: float = Field(ge=0, le=1, description="Surface fraction", default=1.0 / 7)
    emis: float = Field(ge=0, le=1, description="Surface emissivity", default=0.95)
    chanohm: Optional[float] = Field(default=0.0)
    cpanohm: Optional[float] = Field(default=1200.0)
    kkanohm: Optional[float] = Field(default=0.4)
    ohm_threshsw: Optional[float] = Field(default=0.0)
    ohm_threshwd: Optional[float] = Field(default=0.0)
    ohm_coef: Optional[OHMCoefficients] = None
    soildepth: float = Field(default=0.15)
    soilstorecap: float = Field(default=150.0)
    statelimit: float = Field(default=10.0)
    wetthresh: float = Field(default=0.5)
    sathydraulicconduct: float = Field(default=0.0001)
    waterdist: Optional[WaterDistribution] = Field(
        default=None, description="Water distribution parameters"
    )
    storedrainprm: StorageDrainParams = Field(
        default_factory=StorageDrainParams, description="Storage and drain parameters"
    )
    snowpacklimit: Optional[float] = Field(default=10.0)
    thermal_layers: ThermalLayers = Field(
        default_factory=ThermalLayers, description="Thermal layers for the surface"
    )
    irrfrac: Optional[float] = Field(default=0.0)
    _surface_type: Optional[SurfaceType] = PrivateAttr(default=None)

    def set_surface_type(self, surface_type: SurfaceType):
        self._surface_type = surface_type
        if self._surface_type == SurfaceType.WATER:
            if self.waterdist is not None:
                raise ValueError("Water surface should not have water distribution")
        else:
            if self.waterdist is None:
                raise ValueError(
                    f"Water distribution required for {self._surface_type.value}"
                )
            self.waterdist.validate_distribution(self._surface_type)

    def get_surface_type(self) -> SurfaceType:
        return self._surface_type

    def get_surface_name(self) -> str:
        return self._surface_type.value

    def get_surface_index(self) -> int:
        dict_surface_type = {
            SurfaceType.PAVED: 0,
            SurfaceType.BLDGS: 1,
            SurfaceType.EVETR: 2,
            SurfaceType.DECTR: 3,
            SurfaceType.GRASS: 4,
            SurfaceType.BSOIL: 5,
            SurfaceType.WATER: 6,
        }
        return dict_surface_type[self._surface_type]

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert surface properties to DataFrame state format.
        This is the base implementation that handles common surface properties."""
        df_state = init_df_state(grid_id)

        # Get surface index
        surf_idx = self.get_surface_index()

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = f"({surf_idx},)"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = np.nan
            df_state.at[grid_id, (col_name, idx_str)] = value

        # Get all properties of this class using introspection
        properties = [
            attr
            for attr in dir(self)
            if not attr.startswith("_") and not callable(getattr(self, attr))
        ]
        # drop 'surface_type' and model-specific properties
        properties = [
            p for p in properties if p != "surface_type" and not p.startswith("model_")
        ]

        # Process each property
        dfs = [df_state]  # List to collect all DataFrames

        for property in properties:
            # Handle nested properties with their own to_df_state methods
            if property in ["waterdist", "storedrainprm", "thermal_layers", "ohm_coef"]:
                nested_obj = getattr(self, property)
                if nested_obj is not None and hasattr(nested_obj, "to_df_state"):
                    nested_df = nested_obj.to_df_state(grid_id, surf_idx)
                    dfs.append(nested_df)
                continue

            try:
                value = getattr(self, property)
                set_df_value(property, value)
            except Exception as e:
                print(f"Warning: Could not set property {property}: {str(e)}")
                continue

        # Merge all DataFrames
        df_final = pd.concat(dfs, axis=1).sort_index(axis=1)
        return df_final


class NonVegetatedSurfaceProperties(SurfaceProperties):
    alb: float = Field(ge=0, le=1, description="Surface albedo", default=0.1)
    emis: float = Field(ge=0, le=1, description="Surface emissivity", default=0.95)
    z0: float = Field(ge=0, description="Roughness length for momentum", default=0.1)
    zdh: float = Field(ge=0, description="Zero-plane displacement height", default=0.05)
    frfossilfuel_heat: float = Field(
        ge=0, le=1, description="Fraction of fossil fuel heat", default=0.0
    )
    frfossilfuel_cool: float = Field(
        ge=0, le=1, description="Fraction of fossil fuel cooling", default=0.0
    )
    waterdist: Optional[WaterDistribution] = None
    storedrainprm: Optional[StorageDrainParams] = None
    thermal_layers: Optional[ThermalLayer] = None
    ohm_coef: Optional[OHMCoefficients] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert non-vegetated surface properties to DataFrame state format."""
        dfs = []

        # Get base properties from parent
        df_base = super().to_df_state(grid_id)
        dfs.append(df_base)

        surf_idx = self.get_surface_index()

        # Create DataFrame for this class's properties
        param_tuples = []
        values = []

        # Add all non-inherited properties that aren't model-specific or nested objects
        for attr in dir(self):
            if (
                not attr.startswith("_")
                and not callable(getattr(self, attr))
                and not attr.startswith("model_")
                and attr
                not in [
                    "surface_type",
                    "waterdist",
                    "storedrainprm",
                    "thermal_layers",
                    "ohm_coef",
                ]
                and attr not in dir(super())
            ):
                value = getattr(self, attr)
                if not isinstance(value, (BaseModel, Enum)):
                    param_tuples.append((attr, (surf_idx,)))
                    values.append(value)

        if param_tuples:  # Only create DataFrame if we have properties to add
            columns = pd.MultiIndex.from_tuples(param_tuples, names=["var", "ind_dim"])
            df = pd.DataFrame(
                index=pd.Index([grid_id], name="grid"),
                columns=columns,
                data=[values],
                dtype=float,
            )
            dfs.append(df)

        # Add nested property DataFrames
        for nested_prop in ["waterdist", "storedrainprm", "thermal_layers", "ohm_coef"]:
            nested_obj = getattr(self, nested_prop)
            if nested_obj is not None and hasattr(nested_obj, "to_df_state"):
                nested_df = nested_obj.to_df_state(grid_id, surf_idx)
                dfs.append(nested_df)

        # Merge all DataFrames
        df_final = pd.concat(dfs, axis=1)
        return df_final


class PavedProperties(NonVegetatedSurfaceProperties):
    surface_type: Literal[SurfaceType.PAVED] = SurfaceType.PAVED

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert paved surface properties to DataFrame state format."""
        dfs = []

        # Get base properties from parent
        df_base = super().to_df_state(grid_id)
        dfs.append(df_base)

        surf_idx = self.get_surface_index()

        # Create DataFrame for this class's properties
        param_tuples = []
        values = []

        # Add all non-inherited properties that aren't model-specific or nested objects
        for attr in dir(self):
            if (
                not attr.startswith("_")
                and not callable(getattr(self, attr))
                and not attr.startswith("model_")
                and attr
                not in [
                    "surface_type",
                    "waterdist",
                    "storedrainprm",
                    "thermal_layers",
                    "ohm_coef",
                ]
                and attr not in dir(super())
            ):
                value = getattr(self, attr)
                if not isinstance(value, (BaseModel, Enum)):
                    param_tuples.append((attr, (surf_idx,)))
                    values.append(value)

        if param_tuples:  # Only create DataFrame if we have properties to add
            columns = pd.MultiIndex.from_tuples(param_tuples, names=["var", "ind_dim"])
            df = pd.DataFrame(
                index=pd.Index([grid_id], name="grid"),
                columns=columns,
                data=[values],
                dtype=float,
            )
            dfs.append(df)

        # Add nested property DataFrames
        for nested_prop in ["waterdist", "storedrainprm", "thermal_layers", "ohm_coef"]:
            nested_obj = getattr(self, nested_prop)
            if nested_obj is not None and hasattr(nested_obj, "to_df_state"):
                nested_df = nested_obj.to_df_state(grid_id, surf_idx)
                dfs.append(nested_df)

        # Merge all DataFrames
        df_final = pd.concat(dfs, axis=1)
        return df_final


class BuildingLayer(BaseModel):
    alb: float = Field(ge=0, le=1, description="Surface albedo", default=0.1)
    emis: float = Field(ge=0, le=1, description="Surface emissivity", default=0.95)
    thermal_layers: ThermalLayers = Field(
        default_factory=ThermalLayers, description="Thermal layers for the surface"
    )
    statelimit: float = Field(default=10.0)
    soilstorecap: float = Field(default=150.0)
    wetthresh: float = Field(default=0.5)
    roof_albedo_dir_mult_fact: Optional[float] = Field(default=0.1)
    wall_specular_frac: Optional[float] = Field(default=0.1)

    def to_df_state(
        self, grid_id: int, layer_idx: int, is_roof: bool = True
    ) -> pd.DataFrame:
        """Convert building layer parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            layer_idx: Layer index (0 or 1 for two layers)
            is_roof: True if this is a roof layer, False if wall layer

        Returns:
            pd.DataFrame: DataFrame containing building layer parameters
        """
        df_state = init_df_state(grid_id)

        # Determine prefix based on layer type
        prefix = "roof" if is_roof else "wall"

        # Add basic parameters
        df_state[(f"{prefix}_alb", f"({layer_idx},)")] = self.alb
        df_state[(f"{prefix}_emis", f"({layer_idx},)")] = self.emis
        df_state[(f"{prefix}_statelimit", f"({layer_idx},)")] = self.statelimit
        df_state[(f"{prefix}_soilstorecap", f"({layer_idx},)")] = self.soilstorecap
        df_state[(f"{prefix}_wetthresh", f"({layer_idx},)")] = self.wetthresh

        # Add layer-specific parameters
        if is_roof and self.roof_albedo_dir_mult_fact is not None:
            df_state[(f"{prefix}_albedo_dir_mult_fact", f"({layer_idx},)")] = (
                self.roof_albedo_dir_mult_fact
            )
        elif not is_roof and self.wall_specular_frac is not None:
            df_state[(f"{prefix}_specular_frac", f"({layer_idx},)")] = (
                self.wall_specular_frac
            )

        # Add thermal layers
        df_thermal = self.thermal_layers.to_df_state(grid_id, layer_idx)
        df_state = pd.concat([df_state, df_thermal], axis=1)

        return df_state


class VerticalLayers(BaseModel):
    nlayer: int
    height: List[float]
    veg_frac: List[float]
    veg_scale: List[float]
    building_frac: List[float]
    building_scale: List[float]
    roofs: List[BuildingLayer]
    walls: List[BuildingLayer]

    @model_validator(mode="after")
    def validate_building(self) -> "VerticalLayers":
        # Validate building heights
        if len(self.height) != self.nlayer + 1:
            raise ValueError(
                f"Number of building heights ({len(self.height)}) must match nlayer+1 = ({self.nlayer+1})"
            )

        # Validate building fractions
        if len(self.building_frac) != self.nlayer:
            raise ValueError(
                f"Number of building fractions ({len(self.building_frac)}) must match nlayer ({self.nlayer})"
            )
        if not math.isclose(sum(self.building_frac), 1.0, rel_tol=1e-9):
            raise ValueError(
                f"Building fractions must sum to 1.0, got {sum(self.building_frac)}"
            )

        # Validate building scales
        if len(self.building_scale) != self.nlayer:
            raise ValueError(
                f"Number of building scales ({len(self.building_scale)}) must match nlayer ({self.nlayer})"
            )

        # Validate number of roof layers matches nlayer
        if len(self.roofs) != self.nlayer:
            raise ValueError(
                f"Number of roof layers ({len(self.roof)}) must match nlayer ({self.nlayer})"
            )

        # Validate number of wall layers matches nlayer
        if len(self.walls) != self.nlayer:
            raise ValueError(
                f"Number of wall layers ({len(self.wall)}) must match nlayer ({self.nlayer})"
            )

        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert vertical layers to DataFrame state format."""
        # Initialize empty DataFrame with grid_id index
        df_state = init_df_state(grid_id)

        # Set number of vertical layers
        df_state[(f"nlayer", "0")] = self.nlayer

        # Set heights for each layer boundary (nlayer + 1 heights needed)
        for i in range(self.nlayer + 1):
            df_state[("height", f"({i},)")] = self.height[i]

        # Set vegetation and building parameters for each layer
        for var in ["veg_frac", "veg_scale", "building_frac", "building_scale"]:
            for i in range(self.nlayer):
                df_state[(f"{var}", f"({i},)")] = getattr(self, var)[i]

        # Convert roof and wall properties to DataFrame format for each layer
        df_roofs = pd.concat(
            [
                self.roofs[i].to_df_state(grid_id, i, is_roof=True)
                for i in range(self.nlayer)
            ],
            axis=1,
        )
        df_walls = pd.concat(
            [
                self.walls[i].to_df_state(grid_id, i, is_roof=False)
                for i in range(self.nlayer)
            ],
            axis=1,
        )

        # Combine all DataFrames
        df_state = pd.concat([df_state, df_roofs, df_walls], axis=1)

        return df_state


class BuildingProperties(NonVegetatedSurfaceProperties):
    surface_type: Literal[SurfaceType.BLDGS] = SurfaceType.BLDGS
    faibldg: float = Field(
        ge=0, default=0.3, description="Frontal area index of buildings"
    )
    bldgh: float = Field(ge=0, default=10.0, description="Building height")

    @model_validator(mode="after")
    def validate_rsl_zd_range(self) -> "BuildingProperties":
        # Existing validation
        sfr_bldg_lower_limit = 0.18
        if self.sfr < sfr_bldg_lower_limit:
            if self.faibldg < 0.25 * (1 - self.sfr):
                error_message = ValueError(
                    "The Frontal Area Index (FAI) is falling below the lower limit of: 0.25 * (1 - PAI), which is likely causing issues regarding negative displacement height (zd) in the RSL.\n"
                    f"\tYou have entered a building FAI of {self.faibldg} and a building PAI of {self.sfr}.\n"
                    "\tFor more details, please refer to: https://github.com/UMEP-dev/SUEWS/issues/302"
                )
                exceptions.append(error_message)
        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert building properties to DataFrame state format."""
        df_state = super().to_df_state(grid_id)
        surf_idx = self.get_surface_index()

        # Helper function to set values in DataFrame
        # def set_df_value(col_name: str, value: float):
        #     idx_str = f"({surf_idx},)"
        #     if (col_name, idx_str) not in df_state.columns:
        #         df_state[(col_name, idx_str)] = np.nan
        #     df_state.loc[grid_id, (col_name, idx_str)] = value

        def set_df_value(df: pd.DataFrame, grid_id: int, col_name: str, idx_str: str, value: float) -> pd.DataFrame:
            """Helper function to safely set values in DataFrame with MultiIndex columns."""
            col = (col_name, idx_str)
            if col not in df.columns:
                df[col] = np.nan
                df = df.sort_index(axis=1)
            df.loc[grid_id, col] = value
            return df

        # Add all non-inherited properties
        for attr in ["faibldg", "bldgh"]:
            df_state = set_df_value(df_state, grid_id, attr, f"({surf_idx},)", getattr(self, attr))

        return df_state


class BaresoilProperties(NonVegetatedSurfaceProperties):
    surface_type: Literal[SurfaceType.BSOIL] = SurfaceType.BSOIL

    # def to_df_state(self, grid_id: int) -> pd.DataFrame:
    #     """Convert bare soil properties to DataFrame state format."""
    #     df_state = super().to_df_state(grid_id)

    #     return df_state


class WaterProperties(NonVegetatedSurfaceProperties):
    surface_type: Literal[SurfaceType.WATER] = SurfaceType.WATER
    flowchange: float = Field(default=0.0)

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert water surface properties to DataFrame state format."""
        df_state = super().to_df_state(grid_id)
        surf_idx = self.get_surface_index()

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = f"({surf_idx},)"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = np.nan
            df_state.loc[grid_id, (col_name, idx_str)] = value

        list_attr = [
            "flowchange",
        ]

        # Add all non-inherited properties
        for attr in list_attr:
            try:
                value = getattr(self, attr)
                if not isinstance(value, (BaseModel, Enum)):
                    set_df_value(attr, value)
            except Exception as e:
                print(f"Warning: Could not set property {attr}: {str(e)}")

        return df_state


class ModelControl(BaseModel):
    tstep: int = Field(description="Time step in seconds")
    forcing_file: str
    output_file: str
    # daylightsaving_method: int
    diagnose: int

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert model control properties to DataFrame state format."""
        df_state = init_df_state(grid_id)

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = "0"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = np.nan
            df_state.at[grid_id, (col_name, idx_str)] = value

        list_attr = ["tstep", "diagnose"]
        for attr in list_attr:
            set_df_value(attr, getattr(self, attr))
        return df_state


class ModelPhysics(BaseModel):
    netradiationmethod: int
    emissionsmethod: int
    storageheatmethod: int
    ohmincqf: int
    roughlenmommethod: int
    roughlenheatmethod: int
    stabilitymethod: int
    smdmethod: int
    waterusemethod: int
    diagmethod: int
    faimethod: int
    localclimatemethod: int
    snowuse: int

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert model physics properties to DataFrame state format."""
        df_state = init_df_state(grid_id)

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = "0"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = np.nan
            df_state.at[grid_id, (col_name, idx_str)] = int(value)

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
        ]

        for attr in list_attr:
            try:
                properties[attr] = int(df.loc[grid_id, (attr, "0")])
            except KeyError:
                raise ValueError(f"Missing attribute '{attr}' in the DataFrame")

        return cls(**properties)


class LUMPSParams(BaseModel):
    raincover: float = Field(ge=0, le=1, default=0.25)
    rainmaxres: float = Field(ge=0, le=1, default=0.25)
    drainrt: float = Field(ge=0, le=1, default=0.25)
    veg_type: int = Field(default=1)

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
            df_state[(attr, "0")] = getattr(self, attr)

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

        return cls(**params)


class SPARTACUSParams(BaseModel):
    air_ext_lw: float
    air_ext_sw: float
    air_ssa_lw: float
    air_ssa_sw: float
    ground_albedo_dir_mult_fact: float
    n_stream_lw_urban: int
    n_stream_sw_urban: int
    n_vegetation_region_urban: int
    sw_dn_direct_frac: float
    use_sw_direct_albedo: float
    veg_contact_fraction_const: float
    veg_fsd_const: float
    veg_ssa_lw: float
    veg_ssa_sw: float

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
            df_state[(param_name, "0")] = value

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

        params = {param: df.loc[grid_id, (param, "0")] for param in spartacus_params}

        return cls(**params)

    # def to_df_state(self, grid_id: int) -> pd.DataFrame:
    #     """Convert SPARTACUS parameters to DataFrame state format.

    #     Args:
    #         grid_id: Grid ID for the DataFrame index

    #     Returns:
    #         pd.DataFrame: DataFrame containing SPARTACUS parameters
    #     """
    #     df_state = init_df_state(grid_id)

    #     # Add all attributes except private ones
    #     for attr in dir(self):
    #         if not attr.startswith("_") and not callable(getattr(self, attr)):
    #             value = getattr(self, attr)
    #             if not isinstance(value, (BaseModel, Enum)):
    #                 df_state[(attr, "0")] = value

    #     return df_state


class DayProfile(BaseModel):
    working_day: float
    holiday: float

    def to_df_state(self, grid_id: int, param_name: str) -> pd.DataFrame:
        """
        Convert day profile to DataFrame state format.

        Args:
            grid_id (int): Grid ID for the DataFrame index.
            param_name (str): Name of the parameter this profile belongs to.

        Returns:
            pd.DataFrame: DataFrame containing day profile parameters.
        """

        df_state = init_df_state(grid_id)

        day_map = {
            "working_day": 0,
            "holiday": 1,
        }

        for day, idx in day_map.items():
            df_state.loc[grid_id, (param_name, f"({idx},)")] = getattr(self, day)

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, param_name: str
    ) -> "DayProfile":
        """
        Reconstruct DayProfile from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing day profile parameters.
            grid_id (int): Grid ID for the DataFrame index.
            param_name (str): Name of the parameter this profile belongs to.

        Returns:
            DayProfile: Instance of DayProfile.
        """

        day_map = {
            "working_day": 0,
            "holiday": 1,
        }

        # Extract values for working day and holiday from the DataFrame
        params = {}
        for day, idx in day_map.items():
            col = (param_name, f"({idx},)")
            if col in df.columns:
                params[day] = df.loc[grid_id, col]
            else:
                raise KeyError(f"Column {col} not found in DataFrame")

        return cls(**params)

    # # this need to be fixed!
    # def to_df_state(self, grid_id: int, param_name: str) -> pd.DataFrame:
    #     """Convert day profile to DataFrame state format.

    #     Args:
    #         grid_id: Grid ID for the DataFrame index
    #         param_name: Name of the parameter this profile belongs to

    #     Returns:
    #         pd.DataFrame: DataFrame containing day profile parameters
    #     """
    #     df_state = init_df_state(grid_id)

    #     # Set values for working day and holiday
    #     df_state[(f"{param_name}_wd", "0")] = self.working_day
    #     df_state[(f"{param_name}_we", "0")] = self.holiday

    #     return df_state


class WeeklyProfile(BaseModel):
    monday: float
    tuesday: float
    wednesday: float
    thursday: float
    friday: float
    saturday: float
    sunday: float

    def to_df_state(self, grid_id: int, param_name: str) -> pd.DataFrame:
        """Convert weekly profile to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            param_name: Name of the parameter this profile belongs to

        Returns:
            pd.DataFrame: DataFrame containing weekly profile parameters
        """
        df_state = init_df_state(grid_id)

        # Map days to their index
        day_map = {
            "monday": 0,
            "tuesday": 1,
            "wednesday": 2,
            "thursday": 3,
            "friday": 4,
            "saturday": 5,
            "sunday": 6,
        }

        for day, idx in day_map.items():
            df_state[(param_name, f"({idx},)")] = getattr(self, day)

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, param_name: str
    ) -> "WeeklyProfile":
        """Reconstruct WeeklyProfile from a DataFrame state format.

        Args:
            df: DataFrame containing weekly profile parameters
            grid_id: Grid ID for the DataFrame index
            param_name: Name of the parameter to extract values from

        Returns:
            WeeklyProfile: Instance of WeeklyProfile
        """
        # Map days to their index
        day_map = {
            "monday": 0,
            "tuesday": 1,
            "wednesday": 2,
            "thursday": 3,
            "friday": 4,
            "saturday": 5,
            "sunday": 6,
        }

        # Extract values from DataFrame for each day
        params = {
            day: df.loc[grid_id, (param_name, f"({idx},)")]
            for day, idx in day_map.items()
        }

        # Create an instance of WeeklyProfile
        return cls(**params)


class HourlyProfile(BaseModel):
    working_day: Dict[str, float]
    holiday: Dict[str, float]

    @field_validator("working_day", "holiday", mode="before")
    def convert_keys_to_str(cls, v: Dict) -> Dict[str, float]:
        if isinstance(v, dict):
            return {str(k): float(v) for k, v in v.items()}
        return v

    @model_validator(mode="after")
    def validate_hours(self) -> "HourlyProfile":
        for profile in [self.working_day, self.holiday]:
            hours = [int(h) for h in profile.keys()]
            if not all(1 <= h <= 24 for h in hours):
                error_message = ValueError("Hour values must be between 1 and 24")
                exceptions.append(error_message)
                # raise ValueError("Hour values must be between 1 and 24")
            if sorted(hours) != list(range(1, 25)):
                error_message = ValueError("Must have all hours from 1 to 24")
                exceptions.append(error_message)
                # raise ValueError("Must have all hours from 1 to 24")
        return self

    def to_df_state(self, grid_id: int, param_name: str) -> pd.DataFrame:
        """Convert hourly profile to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            param_name: Name of the parameter this profile belongs to

        Returns:
            pd.DataFrame: DataFrame containing hourly profile parameters
        """
        df_state = init_df_state(grid_id)

        # Set working day values (index 0)
        for hour, value in self.working_day.items():
            df_state[(param_name, f"(0,{int(hour)-1})")] = value

        # Set holiday/weekend values (index 1)
        for hour, value in self.holiday.items():
            df_state[(param_name, f"(1,{int(hour)-1})")] = value

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, param_name: str
    ) -> "HourlyProfile":
        """Reconstruct HourlyProfile from a DataFrame state format.

        Args:
            df: DataFrame containing hourly profile parameters
            grid_id: Grid ID for the DataFrame index
            param_name: Name of the parameter to extract values from

        Returns:
            HourlyProfile: Instance of HourlyProfile
        """
        # Extract working day values (index 0)
        working_day = {
            str(hour + 1): df.loc[grid_id, (param_name, f"(0,{hour})")]
            for hour in range(24)
        }

        # Extract holiday/weekend values (index 1)
        holiday = {
            str(hour + 1): df.loc[grid_id, (param_name, f"(1,{hour})")]
            for hour in range(24)
        }

        # Create an instance of HourlyProfile
        return cls(working_day=working_day, holiday=holiday)


class IrrigationParams(BaseModel):
    h_maintain: float
    faut: float
    ie_start: float
    ie_end: float
    internalwateruse_h: float
    daywatper: WeeklyProfile
    daywat: WeeklyProfile
    wuprofa_24hr: HourlyProfile
    wuprofm_24hr: HourlyProfile

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """
        Convert irrigation parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing irrigation parameters
        """

        df_state = init_df_state(grid_id)

        df_state.loc[grid_id, ("h_maintain", 0)] = self.h_maintain
        df_state.loc[grid_id, ("faut", 0)] = self.faut
        df_state.loc[grid_id, ("ie_start", 0)] = self.ie_start
        df_state.loc[grid_id, ("ie_end", 0)] = self.ie_end
        df_state.loc[grid_id, ("internalwateruse_h", 0)] = self.internalwateruse_h

        df_daywatper = self.daywatper.to_df_state(grid_id, "daywatper")
        df_daywat = self.daywat.to_df_state(grid_id, "daywat")

        df_state = df_state.combine_first(df_daywatper)
        df_state = df_state.combine_first(df_daywat)

        df_wuprofa_24hr = self.wuprofa_24hr.to_df_state(grid_id, "wuprofa_24hr")
        df_wuprofm_24hr = self.wuprofm_24hr.to_df_state(grid_id, "wuprofm_24hr")

        df_state = df_state.combine_first(df_wuprofa_24hr)
        df_state = df_state.combine_first(df_wuprofm_24hr)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "IrrigationParams":
        """
        Reconstruct IrrigationParams from a DataFrame state format.

        Args:
            df: DataFrame containing irrigation parameters
            grid_id: Grid ID for the DataFrame index

        Returns:
            IrrigationParams: Instance of IrrigationParams
        """
        # Extract scalar attributes
        h_maintain = df.loc[grid_id, ("h_maintain", 0)]
        faut = df.loc[grid_id, ("faut", 0)]
        ie_start = df.loc[grid_id, ("ie_start", 0)]
        ie_end = df.loc[grid_id, ("ie_end", 0)]
        internalwateruse_h = df.loc[grid_id, ("internalwateruse_h", 0)]

        # Extract WeeklyProfile attributes
        daywatper = WeeklyProfile.from_df_state(df, grid_id, "daywatper")
        daywat = WeeklyProfile.from_df_state(df, grid_id, "daywat")

        # Extract HourlyProfile attributes
        wuprofa_24hr = HourlyProfile.from_df_state(df, grid_id, "wuprofa_24hr")
        wuprofm_24hr = HourlyProfile.from_df_state(df, grid_id, "wuprofm_24hr")

        # Construct and return the IrrigationParams instance
        return cls(
            h_maintain=h_maintain,
            faut=faut,
            ie_start=ie_start,
            ie_end=ie_end,
            internalwateruse_h=internalwateruse_h,
            daywatper=daywatper,
            daywat=daywat,
            wuprofa_24hr=wuprofa_24hr,
            wuprofm_24hr=wuprofm_24hr,
        )


class AnthropogenicHeat(BaseModel):
    qf0_beu: DayProfile
    qf_a: DayProfile
    qf_b: DayProfile
    qf_c: DayProfile
    baset_cooling: DayProfile
    baset_heating: DayProfile
    ah_min: DayProfile
    ah_slope_cooling: DayProfile
    ah_slope_heating: DayProfile
    ahprof_24hr: HourlyProfile
    popdensdaytime: DayProfile
    popdensnighttime: float
    popprof_24hr: HourlyProfile

    # DayProfile coulmns need to be fixed
    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """
        Convert anthropogenic heat parameters to DataFrame state format.

        Args:
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            pd.DataFrame: DataFrame containing anthropogenic heat parameters.
        """

        df_state = init_df_state(grid_id)

        day_profiles = {
            "qf0_beu": self.qf0_beu,
            "qf_a": self.qf_a,
            "qf_b": self.qf_b,
            "qf_c": self.qf_c,
            "baset_cooling": self.baset_cooling,
            "baset_heating": self.baset_heating,
            "ah_min": self.ah_min,
            "ah_slope_cooling": self.ah_slope_cooling,
            "ah_slope_heating": self.ah_slope_heating,
            "popdensdaytime": self.popdensdaytime,
        }
        for param_name, profile in day_profiles.items():
            df_day_profile = profile.to_df_state(grid_id, param_name)
            df_state = df_state.combine_first(df_day_profile)

        hourly_profiles = {
            "ahprof_24hr": self.ahprof_24hr,
            "popprof_24hr": self.popprof_24hr,
        }
        for param_name, profile in hourly_profiles.items():
            df_hourly_profile = profile.to_df_state(grid_id, param_name)
            df_state = df_state.combine_first(df_hourly_profile)

        df_state.loc[grid_id, ("popdensnighttime", 0)] = self.popdensnighttime

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "AnthropogenicHeat":
        """
        Reconstruct AnthropogenicHeat from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing anthropogenic heat parameters.
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            AnthropogenicHeat: Instance of AnthropogenicHeat.
        """

        # Extract DayProfile attributes
        day_profiles = {
            "qf0_beu": DayProfile.from_df_state(df, grid_id, "qf0_beu"),
            "qf_a": DayProfile.from_df_state(df, grid_id, "qf_a"),
            "qf_b": DayProfile.from_df_state(df, grid_id, "qf_b"),
            "qf_c": DayProfile.from_df_state(df, grid_id, "qf_c"),
            "baset_cooling": DayProfile.from_df_state(df, grid_id, "baset_cooling"),
            "baset_heating": DayProfile.from_df_state(df, grid_id, "baset_heating"),
            "ah_min": DayProfile.from_df_state(df, grid_id, "ah_min"),
            "ah_slope_cooling": DayProfile.from_df_state(
                df, grid_id, "ah_slope_cooling"
            ),
            "ah_slope_heating": DayProfile.from_df_state(
                df, grid_id, "ah_slope_heating"
            ),
            "popdensdaytime": DayProfile.from_df_state(df, grid_id, "popdensdaytime"),
        }

        # Extract HourlyProfile attributes
        hourly_profiles = {
            "ahprof_24hr": HourlyProfile.from_df_state(df, grid_id, "ahprof_24hr"),
            "popprof_24hr": HourlyProfile.from_df_state(df, grid_id, "popprof_24hr"),
        }

        # Extract scalar attribute
        popdensnighttime = df.loc[grid_id, ("popdensnighttime", 0)]

        # Construct and return AnthropogenicHeat instance
        return cls(
            **day_profiles,
            **hourly_profiles,
            popdensnighttime=popdensnighttime,
        )


class CO2Params(BaseModel):
    co2pointsource: float
    ef_umolco2perj: float
    enef_v_jkm: float
    fcef_v_kgkm: DayProfile
    frfossilfuel_heat: float
    frfossilfuel_nonheat: float
    maxfcmetab: float
    maxqfmetab: float
    minfcmetab: float
    minqfmetab: float
    trafficrate: DayProfile
    trafficunits: float
    traffprof_24hr: HourlyProfile
    humactivity_24hr: HourlyProfile

    # DayProfile coulmns need to be fixed
    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """
        Convert CO2 parameters to DataFrame state format.

        Args:
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            pd.DataFrame: DataFrame containing CO2 parameters.
        """

        df_state = init_df_state(grid_id)

        scalar_params = {
            "co2pointsource": self.co2pointsource,
            "ef_umolco2perj": self.ef_umolco2perj,
            "enef_v_jkm": self.enef_v_jkm,
            "frfossilfuel_heat": self.frfossilfuel_heat,
            "frfossilfuel_nonheat": self.frfossilfuel_nonheat,
            "maxfcmetab": self.maxfcmetab,
            "maxqfmetab": self.maxqfmetab,
            "minfcmetab": self.minfcmetab,
            "minqfmetab": self.minqfmetab,
            "trafficunits": self.trafficunits,
        }
        for param_name, value in scalar_params.items():
            df_state.loc[grid_id, (param_name, 0)] = value

        day_profiles = {
            "fcef_v_kgkm": self.fcef_v_kgkm,
            "trafficrate": self.trafficrate,
        }
        for param_name, profile in day_profiles.items():
            df_day_profile = profile.to_df_state(grid_id, param_name)
            df_state = df_state.combine_first(df_day_profile)

        hourly_profiles = {
            "traffprof_24hr": self.traffprof_24hr,
            "humactivity_24hr": self.humactivity_24hr,
        }
        for param_name, profile in hourly_profiles.items():
            df_hourly_profile = profile.to_df_state(grid_id, param_name)
            df_state = df_state.combine_first(df_hourly_profile)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "CO2Params":
        """
        Reconstruct CO2Params from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing CO2 parameters.
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            CO2Params: Instance of CO2Params.
        """

        # Extract scalar attributes
        scalar_params = {
            "co2pointsource": df.loc[grid_id, ("co2pointsource", 0)],
            "ef_umolco2perj": df.loc[grid_id, ("ef_umolco2perj", 0)],
            "enef_v_jkm": df.loc[grid_id, ("enef_v_jkm", 0)],
            "frfossilfuel_heat": df.loc[grid_id, ("frfossilfuel_heat", 0)],
            "frfossilfuel_nonheat": df.loc[grid_id, ("frfossilfuel_nonheat", 0)],
            "maxfcmetab": df.loc[grid_id, ("maxfcmetab", 0)],
            "maxqfmetab": df.loc[grid_id, ("maxqfmetab", 0)],
            "minfcmetab": df.loc[grid_id, ("minfcmetab", 0)],
            "minqfmetab": df.loc[grid_id, ("minqfmetab", 0)],
            "trafficunits": df.loc[grid_id, ("trafficunits", 0)],
        }

        # Extract DayProfile attributes
        day_profiles = {
            "fcef_v_kgkm": DayProfile.from_df_state(df, grid_id, "fcef_v_kgkm"),
            "trafficrate": DayProfile.from_df_state(df, grid_id, "trafficrate"),
        }

        # Extract HourlyProfile attributes
        hourly_profiles = {
            "traffprof_24hr": HourlyProfile.from_df_state(
                df, grid_id, "traffprof_24hr"
            ),
            "humactivity_24hr": HourlyProfile.from_df_state(
                df, grid_id, "humactivity_24hr"
            ),
        }

        # Construct and return CO2Params instance
        return cls(
            **scalar_params,
            **day_profiles,
            **hourly_profiles,
        )


class AnthropogenicEmissions(BaseModel):
    startdls: float
    enddls: float
    heat: AnthropogenicHeat
    co2: CO2Params

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """
        Convert anthropogenic emissions parameters to DataFrame state format.

        Args:
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            pd.DataFrame: DataFrame containing anthropogenic emissions parameters.
        """
        df_state = init_df_state(grid_id)

        # Set start and end daylight saving times
        df_state.loc[grid_id, ("startdls", 0)] = self.startdls
        df_state.loc[grid_id, ("enddls", 0)] = self.enddls

        # Add heat parameters
        df_heat = self.heat.to_df_state(grid_id)
        df_state = pd.concat([df_state, df_heat], axis=1)

        # Add CO2 parameters
        df_co2 = self.co2.to_df_state(grid_id)
        df_state = pd.concat([df_state, df_co2], axis=1)

        # Drop duplicate columns if necessary
        df_state = df_state.loc[:, ~df_state.columns.duplicated()]

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "AnthropogenicEmissions":
        """
        Reconstruct AnthropogenicEmissions from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing anthropogenic emissions parameters.
            grid_id (int): Grid ID for the DataFrame index.

        Returns:
            AnthropogenicEmissions: Instance of AnthropogenicEmissions.
        """
        startdls = df.loc[grid_id, ("startdls", 0)]
        enddls = df.loc[grid_id, ("enddls", 0)]

        # Reconstruct heat parameters
        heat = AnthropogenicHeat.from_df_state(df, grid_id)

        # Reconstruct CO2 parameters
        co2 = CO2Params.from_df_state(df, grid_id)

        return cls(startdls=startdls, enddls=enddls, heat=heat, co2=co2)


class Conductance(BaseModel):
    g_max: float = Field(description="Maximum conductance")
    g_k: float = Field(
        description="Conductance parameter related to incoming solar radiation"
    )
    g_q_base: float = Field(
        description="Base value for conductance parameter related to vapor pressure deficit"
    )
    g_q_shape: float = Field(
        description="Shape parameter for conductance related to vapor pressure deficit"
    )
    g_t: float = Field(description="Conductance parameter related to air temperature")
    g_sm: float = Field(description="Conductance parameter related to soil moisture")
    kmax: float = Field(description="Maximum incoming shortwave radiation")
    gsmodel: int = Field(description="Stomatal conductance model selection")
    s1: float = Field(description="Soil moisture threshold parameter")
    s2: float = Field(description="Soil moisture threshold parameter")
    tl: float = Field(description="Air temperature threshold parameter")
    th: float = Field(description="Air temperature threshold parameter")

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
            "gsmodel": self.gsmodel,
            "s1": self.s1,
            "s2": self.s2,
            "tl": self.tl,
            "th": self.th,
        }

        for param_name, value in scalar_params.items():
            df_state.loc[grid_id, (param_name, 0)] = value

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
            "g_max": df.loc[grid_id, ("g_max", 0)],
            "g_k": df.loc[grid_id, ("g_k", 0)],
            "g_q_base": df.loc[grid_id, ("g_q_base", 0)],
            "g_q_shape": df.loc[grid_id, ("g_q_shape", 0)],
            "g_t": df.loc[grid_id, ("g_t", 0)],
            "g_sm": df.loc[grid_id, ("g_sm", 0)],
            "kmax": df.loc[grid_id, ("kmax", 0)],
            "gsmodel": int(df.loc[grid_id, ("gsmodel", 0)]),
            "s1": df.loc[grid_id, ("s1", 0)],
            "s2": df.loc[grid_id, ("s2", 0)],
            "tl": df.loc[grid_id, ("tl", 0)],
            "th": df.loc[grid_id, ("th", 0)],
        }

        return cls(**scalar_params)


class LAIPowerCoefficients(BaseModel):
    growth_lai: float = Field(
        default=0.1,
        description="Power coefficient for LAI in growth equation (LAIPower[1])",
    )
    growth_gdd: float = Field(
        default=0.1,
        description="Power coefficient for GDD in growth equation (LAIPower[2])",
    )
    senescence_lai: float = Field(
        default=0.1,
        description="Power coefficient for LAI in senescence equation (LAIPower[3])",
    )
    senescence_sdd: float = Field(
        default=0.1,
        description="Power coefficient for SDD in senescence equation (LAIPower[4])",
    )

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
                df_state[(col_name, idx_str)] = np.nan
            df_state.at[grid_id, (col_name, idx_str)] = value

        # Set power coefficients in order
        for i, value in enumerate(self.to_list()):
            set_df_value("laipower", (i, veg_idx), value)

        return df_state


class LAIParams(BaseModel):
    baset: float = Field(
        default=10.0,
        description="Base Temperature for initiating growing degree days (GDD) for leaf growth [degC]",
    )
    gddfull: float = Field(
        default=100.0,
        description="Growing degree days (GDD) needed for full capacity of LAI [degC]",
    )
    basete: float = Field(
        default=10.0,
        description="Base temperature for initiating senescence degree days (SDD) for leaf off [degC]",
    )
    sddfull: float = Field(
        default=100.0,
        description="Senescence degree days (SDD) needed to initiate leaf off [degC]",
    )
    laimin: float = Field(default=0.1, description="Leaf-off wintertime value [m2 m-2]")
    laimax: float = Field(
        default=10.0, description="Full leaf-on summertime value [m2 m-2]"
    )
    laipower: LAIPowerCoefficients = Field(
        default_factory=LAIPowerCoefficients,
        description="LAI calculation power parameters for growth and senescence",
    )
    laitype: int = Field(
        default=0,
        description="LAI calculation choice (0: original, 1: new high latitude)",
    )

    @model_validator(mode="after")
    def validate_lai_ranges(self) -> "LAIParams":
        if self.laimin > self.laimax:
            error_message = ValueError(
                f"laimin ({self.laimin}) must be less than or equal to laimax ({self.laimax})."
            )
            exceptions.append(error_message)
            # raise ValueError(f"laimin ({self.laimin})must be less than or equal to laimax ({self.laimax}).")
        if self.baset > self.gddfull:
            error_message = ValueError(
                f"baset ({self.baset}) must be less than gddfull ({self.gddfull})."
            )
            exceptions.append(error_message)
            # raise ValueError(f"baset {self.baset} must be less than gddfull ({self.gddfull}).")
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
                df_state[(col_name, idx_str)] = np.nan
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
            set_df_value(param, (veg_idx,), value)

        # Add LAI power coefficients using the LAIPowerCoefficients to_df_state method
        if self.laipower:
            df_power = self.laipower.to_df_state(grid_id, veg_idx)
            # Merge power coefficients into main DataFrame
            for col in df_power.columns:
                if col[0] != "grid_iv":  # Skip the grid_iv column
                    df_state[col] = df_power[col]

        return df_state


class VegetatedSurfaceProperties(SurfaceProperties):
    alb_min: float = Field(ge=0, le=1, description="Minimum albedo", default=0.2)
    alb_max: float = Field(ge=0, le=1, description="Maximum albedo", default=0.3)
    beta_bioco2: float = Field(
        default=0.6, description="Biogenic CO2 exchange coefficient"
    )
    beta_enh_bioco2: float = Field(
        default=0.7, description="Enhanced biogenic CO2 exchange coefficient"
    )
    alpha_bioco2: float = Field(
        default=0.8, description="Biogenic CO2 exchange coefficient"
    )
    alpha_enh_bioco2: float = Field(
        default=0.9, description="Enhanced biogenic CO2 exchange coefficient"
    )
    resp_a: float = Field(default=1.0, description="Respiration coefficient")
    resp_b: float = Field(default=1.1, description="Respiration coefficient")
    theta_bioco2: float = Field(
        default=1.2, description="Biogenic CO2 exchange coefficient"
    )
    maxconductance: float = Field(
        default=0.5, description="Maximum surface conductance"
    )
    min_res_bioco2: float = Field(
        default=0.1, description="Minimum respiratory biogenic CO2"
    )
    lai: LAIParams = Field(
        default_factory=LAIParams, description="Leaf area index parameters"
    )
    ie_a: float = Field(
        default=0.5, description="Irrigation efficiency coefficient-automatic"
    )
    ie_m: float = Field(
        default=0.6, description="Irrigation efficiency coefficient-manual"
    )

    @model_validator(mode="after")
    def validate_albedo_range(self) -> "VegetatedSurfaceProperties":
        if self.alb_min > self.alb_max:
            error_message = ValueError(
                f"alb_min (input {self.alb_min}) must be less than or equal to alb_max (entered {self.alb_max})."
            )
            exceptions.append(error_message)
            # raise ValueError(f"alb_min (input {self.alb_min}) must be less than or equal to alb_max (entered {self.alb_max}).")
        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert vegetated surface properties to DataFrame state format."""
        # Get base properties
        df_state = super().to_df_state(grid_id)

        # Add vegetation-specific properties
        surf_idx = self.get_surface_index()

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, idx_str: str, value: float):
            print(f"Setting {col_name} with index {idx_str} to {value}")
            if (col_name, idx_str) not in df_state.columns:
                print(f"Column {col_name} with index {idx_str} not found, adding nan")
                df_state[(col_name, idx_str)] = np.nan
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # Add LAI-related properties
        if hasattr(self, "lai"):
            lai = self.lai
            idx = surf_idx - 2  # Adjust index for vegetation surfaces
            for i, var in enumerate(
                [
                    "albmax",
                    "albmin",
                    "lai_max",
                    "lai_min",
                    "baset",
                    "basete",
                    "gdd_id",
                    "sddfull",
                    "senescencesdd",
                ]
            ):
                if hasattr(lai.laipower, var):
                    set_df_value("laipower", f"{(i, idx)}", getattr(lai.laipower, var))

        # Add CO2-related properties for vegetated surfaces
        if hasattr(self, "beta_bioco2"):
            set_df_value("beta_bioco2", f"({idx},)", self.beta_bioco2)
            set_df_value("beta_enh_bioco2", f"({idx},)", self.beta_enh_bioco2)
            set_df_value("alpha_bioco2", f"({idx},)", self.alpha_bioco2)
            set_df_value("alpha_enh_bioco2", f"({idx},)", self.alpha_enh_bioco2)

        return df_state


class DectrProperties(VegetatedSurfaceProperties):
    faidectree: float = Field(
        default=0.1, description="Frontal area index of deciduous trees"
    )
    dectreeh: float = Field(default=15.0, description="Deciduous tree height")
    pormin_dec: float = Field(default=0.2, description="Minimum porosity")
    pormax_dec: float = Field(default=0.6, description="Maximum porosity")

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert deciduous tree properties to DataFrame state format."""
        # Get base properties from parent
        df_state = super().to_df_state(grid_id)
        surf_idx = self.get_surface_index()

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = f"({surf_idx},)"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = np.nan
            df_state.at[grid_id, (col_name, idx_str)] = value

        # Add all non-inherited properties
        for attr in dir(self):
            if (
                not attr.startswith("_")
                and not callable(getattr(self, attr))
                and attr not in ["surface_type"]
                and attr not in dir(super())
            ):
                try:
                    value = getattr(self, attr)
                    if not isinstance(value, (BaseModel, Enum)):
                        set_df_value(attr, value)
                except Exception as e:
                    print(f"Warning: Could not set property {attr}: {str(e)}")

        return df_state


class EvetrProperties(VegetatedSurfaceProperties):
    faievetree: float = Field(
        default=0.1, description="Frontal area index of evergreen trees"
    )
    evetreeh: float = Field(default=15.0, description="Evergreen tree height")

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert evergreen tree properties to DataFrame state format."""
        # Get base properties from parent
        df_state = super().to_df_state(grid_id)
        surf_idx = self.get_surface_index()

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = f"({surf_idx},)"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = np.nan
            df_state.at[grid_id, (col_name, idx_str)] = value

        # Add all non-inherited properties
        for attr in dir(self):
            if (
                not attr.startswith("_")
                and not callable(getattr(self, attr))
                and attr not in ["surface_type"]
                and attr not in dir(super())
            ):
                try:
                    value = getattr(self, attr)
                    if not isinstance(value, (BaseModel, Enum)):
                        set_df_value(attr, value)
                except Exception as e:
                    print(f"Warning: Could not set property {attr}: {str(e)}")

        return df_state


class SnowParams(BaseModel):
    crwmax: float
    crwmin: float
    narp_emis_snow: float
    preciplimit: float
    preciplimitalb: float
    snowalbmax: float
    snowalbmin: float
    snowdensmin: float
    snowdensmax: float
    snowlimbldg: float
    snowlimpaved: float
    snowprof_24hr: HourlyProfile
    tau_a: float
    tau_f: float
    tau_r: float
    tempmeltfact: float
    radmeltfact: float

    @model_validator(mode="after")
    def validate_crw_range(self) -> "SnowParams":
        if self.crwmin >= self.crwmax:
            error_message = ValueError(
                f"crwmin ({self.crwmin}) must be less than crwmax ({self.crwmax})."
            )
            exceptions.append(error_message)
            # raise ValueError(f"crwmin ({self.crwmin}) must be less than crwmax ({self.crwmax}).")
        return self

    @model_validator(mode="after")
    def validate_snowalb_range(self) -> "SnowParams":
        if self.snowalbmin >= self.snowalbmax:
            error_message = ValueError(
                f"snowalbmin ({self.snowalbmin}) must be less than snowalbmax ({self.snowalbmax})."
            )
            exceptions.append(error_message)
            # raise ValueError(f"snowalbmin ({self.snowalbmin}) must be less than snowalbmax ({self.snowalbmax}).")
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
            df_state.loc[grid_id, (param_name, 0)] = value

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
            "crwmax": df.loc[grid_id, ("crwmax", 0)],
            "crwmin": df.loc[grid_id, ("crwmin", 0)],
            "narp_emis_snow": df.loc[grid_id, ("narp_emis_snow", 0)],
            "preciplimit": df.loc[grid_id, ("preciplimit", 0)],
            "preciplimitalb": df.loc[grid_id, ("preciplimitalb", 0)],
            "snowalbmax": df.loc[grid_id, ("snowalbmax", 0)],
            "snowalbmin": df.loc[grid_id, ("snowalbmin", 0)],
            "snowdensmin": df.loc[grid_id, ("snowdensmin", 0)],
            "snowdensmax": df.loc[grid_id, ("snowdensmax", 0)],
            "snowlimbldg": df.loc[grid_id, ("snowlimbldg", 0)],
            "snowlimpaved": df.loc[grid_id, ("snowlimpaved", 0)],
            "tau_a": df.loc[grid_id, ("tau_a", 0)],
            "tau_f": df.loc[grid_id, ("tau_f", 0)],
            "tau_r": df.loc[grid_id, ("tau_r", 0)],
            "tempmeltfact": df.loc[grid_id, ("tempmeltfact", 0)],
            "radmeltfact": df.loc[grid_id, ("radmeltfact", 0)],
        }

        # Extract HourlyProfile
        snowprof_24hr = HourlyProfile.from_df_state(df, grid_id, "snowprof_24hr")

        # Construct and return the SnowParams instance
        return cls(snowprof_24hr=snowprof_24hr, **scalar_params)


class LandCover(BaseModel):
    paved: PavedProperties
    bldgs: BuildingProperties
    dectr: DectrProperties
    evetr: EvetrProperties
    grass: VegetatedSurfaceProperties
    bsoil: BaresoilProperties
    water: WaterProperties

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

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert land cover to DataFrame state format"""
        # df_state = init_df_state(grid_id)

        list_df_state = []
        for lc in ["paved", "bldgs", "dectr", "evetr", "grass", "bsoil", "water"]:
            print('')
            print(lc)
            df_state = getattr(self, lc).to_df_state(grid_id)
            list_df_state.append(df_state)
        df_state = pd.concat(list_df_state, axis=1)
        return df_state


class SiteProperties(BaseModel):
    lat: float = Field(ge=-90, le=90)
    lng: float = Field(ge=-180, le=180)
    alt: float
    timezone: int = Field(ge=-12, le=12)
    surfacearea: float = Field(gt=0)
    z: float = Field(gt=0)
    z0m_in: float = Field(gt=0)
    zdm_in: float = Field(gt=0)
    pipecapacity: float = Field(gt=0)
    runofftowater: float = Field(ge=0, le=1)
    narp_trans_site: float
    lumps: LUMPSParams
    spartacus: SPARTACUSParams
    conductance: Conductance
    irrigation: IrrigationParams
    anthropogenic_emissions: AnthropogenicEmissions
    snow: SnowParams
    land_cover: LandCover
    vertical_layers: VerticalLayers

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert site properties to DataFrame state format"""
        df_state = init_df_state(grid_id)

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
        ]:
            df_state.loc[grid_id, f"({var}, 0)"] = getattr(self, var)


        return df_state


class Site(BaseModel):
    name: str
    gridiv: int
    properties: SiteProperties
    initial_states: InitialStates

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert site to DataFrame state format"""
        df_site = self.properties.to_df_state(grid_id)
        df_initial_states = self.initial_states.to_df_state(grid_id)
        df_state = pd.concat([df_site, df_initial_states], axis=1)
        return df_state


class Model(BaseModel):
    control: ModelControl
    physics: ModelPhysics


class SUEWSConfig(BaseModel):
    name: str
    description: str
    model: Model
    site: List[Site]

    class Config:
        extra = "allow"

    def create_multi_index_columns(self, columns_file: str) -> pd.MultiIndex:
        """Create MultiIndex from df_state_columns.txt"""
        with open(columns_file, "r") as f:
            lines = f.readlines()

        tuples = []
        for line in lines:
            col_name, indices = line.strip().split(",", 1)
            str_indices = f"{indices}" if indices != "0" else "0"
            tuples.append((col_name, str_indices))

        return pd.MultiIndex.from_tuples(tuples)

    def to_df_state(self) -> pd.DataFrame:
        """Convert config to DataFrame state format"""
        # Initialize empty DataFrame with correct structure
        columns = self.create_multi_index_columns("df_state_columns.txt")
        columns.set_names(["var", "ind_dim"], inplace=True)
        index = pd.Index([0], name="grid")
        df = pd.DataFrame(index=index, columns=columns)

        # Process each site (assuming single site for now)
        site = self.site[0]
        props = site.properties
        gridiv = site.gridiv
        # update index using gridiv
        df.index = pd.Index([gridiv], name="grid")

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, indices: Union[int, Tuple], value: float):
            if isinstance(indices, int):
                if indices == 0:
                    str_indices = str(indices)
                else:
                    str_indices = f"({indices},)"
            else:
                # Tuples should maintain their string representation
                str_indices = str(indices)
            df.loc[gridiv, (col_name, str_indices)] = value

        # Model control
        set_df_value("tstep", 0, self.model.control.tstep)
        set_df_value("diagnose", 0, self.model.control.diagnose)
        # Model physics
        set_df_value("netradiationmethod", 0, self.model.physics.netradiationmethod)
        set_df_value("emissionsmethod", 0, self.model.physics.emissionsmethod)
        set_df_value("storageheatmethod", 0, self.model.physics.storageheatmethod)
        set_df_value("ohmincqf", 0, self.model.physics.ohmincqf)
        set_df_value("roughlenmommethod", 0, self.model.physics.roughlenmommethod)
        set_df_value("roughlenheatmethod", 0, self.model.physics.roughlenheatmethod)
        set_df_value("stabilitymethod", 0, self.model.physics.stabilitymethod)
        set_df_value("smdmethod", 0, self.model.physics.smdmethod)
        set_df_value("waterusemethod", 0, self.model.physics.waterusemethod)
        set_df_value("diagmethod", 0, self.model.physics.diagmethod)
        set_df_value("faimethod", 0, self.model.physics.faimethod)
        set_df_value("localclimatemethod", 0, self.model.physics.localclimatemethod)
        set_df_value("snowuse", 0, self.model.physics.snowuse)

        # Basic site properties
        set_df_value("gridiv", 0, gridiv)
        set_df_value("lat", 0, props.lat)
        set_df_value("lng", 0, props.lng)
        set_df_value("alt", 0, props.alt)
        set_df_value("timezone", 0, props.timezone)
        set_df_value("surfacearea", 0, props.surfacearea)
        set_df_value("z", 0, props.z)
        set_df_value("z0m_in", 0, props.z0m_in)
        set_df_value("zdm_in", 0, props.zdm_in)
        set_df_value("pipecapacity", 0, props.pipecapacity)
        set_df_value("runofftowater", 0, props.runofftowater)
        set_df_value("narp_trans_site", 0, props.narp_trans_site)

        # LUMPS parameters
        lumps = props.lumps
        set_df_value("raincover", 0, lumps.raincover)
        set_df_value("rainmaxres", 0, lumps.rainmaxres)
        set_df_value("drainrt", 0, lumps.drainrt)
        set_df_value("veg_type", 0, lumps.veg_type)

        # SPARTACUS parameters
        spartacus = props.spartacus
        set_df_value("air_ext_lw", 0, spartacus.air_ext_lw)
        set_df_value("air_ext_sw", 0, spartacus.air_ext_sw)
        set_df_value("air_ssa_lw", 0, spartacus.air_ssa_lw)
        set_df_value("air_ssa_sw", 0, spartacus.air_ssa_sw)
        set_df_value(
            "ground_albedo_dir_mult_fact", 0, spartacus.ground_albedo_dir_mult_fact
        )
        set_df_value("n_stream_lw_urban", 0, spartacus.n_stream_lw_urban)
        set_df_value("n_stream_sw_urban", 0, spartacus.n_stream_sw_urban)
        set_df_value(
            "n_vegetation_region_urban", 0, spartacus.n_vegetation_region_urban
        )
        set_df_value("sw_dn_direct_frac", 0, spartacus.sw_dn_direct_frac)
        set_df_value("use_sw_direct_albedo", 0, spartacus.use_sw_direct_albedo)
        set_df_value(
            "veg_contact_fraction_const", 0, spartacus.veg_contact_fraction_const
        )
        set_df_value("veg_fsd_const", 0, spartacus.veg_fsd_const)
        set_df_value("veg_ssa_lw", 0, spartacus.veg_ssa_lw)
        set_df_value("veg_ssa_sw", 0, spartacus.veg_ssa_sw)

        # Conductance parameters
        cond = props.conductance
        set_df_value("g_max", 0, cond.g_max)
        set_df_value("g_k", 0, cond.g_k)
        set_df_value("g_q_base", 0, cond.g_q_base)
        set_df_value("g_q_shape", 0, cond.g_q_shape)
        set_df_value("g_t", 0, cond.g_t)
        set_df_value("g_sm", 0, cond.g_sm)
        set_df_value("kmax", 0, cond.kmax)
        set_df_value("gsmodel", 0, cond.gsmodel)
        set_df_value("s1", 0, cond.s1)
        set_df_value("s2", 0, cond.s2)

        # Irrigation parameters
        irr = props.irrigation
        set_df_value("h_maintain", 0, irr.h_maintain)
        set_df_value("faut", 0, irr.faut)
        set_df_value("ie_start", 0, irr.ie_start)
        set_df_value("ie_end", 0, irr.ie_end)
        set_df_value("internalwateruse_h", 0, irr.internalwateruse_h)

        # Daily water parameters
        set_df_value("daywatper", (0,), irr.daywatper.monday)
        set_df_value("daywatper", (1,), irr.daywatper.tuesday)
        set_df_value("daywatper", (2,), irr.daywatper.wednesday)
        set_df_value("daywatper", (3,), irr.daywatper.thursday)
        set_df_value("daywatper", (4,), irr.daywatper.friday)
        set_df_value("daywatper", (5,), irr.daywatper.saturday)
        set_df_value("daywatper", (6,), irr.daywatper.sunday)

        set_df_value("daywat", (0,), irr.daywat.monday)
        set_df_value("daywat", (1,), irr.daywat.tuesday)
        set_df_value("daywat", (2,), irr.daywat.wednesday)
        set_df_value("daywat", (3,), irr.daywat.thursday)
        set_df_value("daywat", (4,), irr.daywat.friday)
        set_df_value("daywat", (5,), irr.daywat.saturday)
        set_df_value("daywat", (6,), irr.daywat.sunday)

        # Water use profile
        for hour in range(24):
            hour_str = str(hour + 1)
            set_df_value(
                "wuprofa_24hr", (hour, 0), irr.wuprofa_24hr.working_day[hour_str]
            )
            set_df_value("wuprofa_24hr", (hour, 1), irr.wuprofa_24hr.holiday[hour_str])

        # SPARTACUS parameters
        spartacus = props.spartacus
        set_df_value("air_ext_lw", 0, spartacus.air_ext_lw)
        set_df_value("air_ext_sw", 0, spartacus.air_ext_sw)
        set_df_value("air_ssa_lw", 0, spartacus.air_ssa_lw)
        set_df_value("air_ssa_sw", 0, spartacus.air_ssa_sw)
        set_df_value(
            "ground_albedo_dir_mult_fact", 0, spartacus.ground_albedo_dir_mult_fact
        )
        set_df_value("n_stream_lw_urban", 0, spartacus.n_stream_lw_urban)
        set_df_value("n_stream_sw_urban", 0, spartacus.n_stream_sw_urban)
        set_df_value(
            "n_vegetation_region_urban", 0, spartacus.n_vegetation_region_urban
        )
        set_df_value("sw_dn_direct_frac", 0, spartacus.sw_dn_direct_frac)
        set_df_value("use_sw_direct_albedo", 0, spartacus.use_sw_direct_albedo)
        set_df_value(
            "veg_contact_fraction_const", 0, spartacus.veg_contact_fraction_const
        )
        set_df_value("veg_fsd_const", 0, spartacus.veg_fsd_const)
        set_df_value("veg_ssa_lw", 0, spartacus.veg_ssa_lw)
        set_df_value("veg_ssa_sw", 0, spartacus.veg_ssa_sw)

        # Conductance properties
        conductance = props.conductance
        set_df_value("g_max", 0, conductance.g_max)
        set_df_value("g_k", 0, conductance.g_k)
        set_df_value("g_q_base", 0, conductance.g_q_base)
        set_df_value("g_q_shape", 0, conductance.g_q_shape)
        set_df_value("g_t", 0, conductance.g_t)
        set_df_value("g_sm", 0, conductance.g_sm)
        set_df_value("kmax", 0, conductance.kmax)
        set_df_value("gsmodel", 0, conductance.gsmodel)
        set_df_value("s1", 0, conductance.s1)
        set_df_value("s2", 0, conductance.s2)
        set_df_value("tl", 0, conductance.tl)
        set_df_value("th", 0, conductance.th)

        # Irrigation parameters
        irrigation = props.irrigation
        set_df_value("h_maintain", 0, irrigation.h_maintain)
        set_df_value("faut", 0, irrigation.faut)
        set_df_value("ie_start", 0, irrigation.ie_start)
        set_df_value("ie_end", 0, irrigation.ie_end)
        set_df_value("internalwateruse_h", 0, irrigation.internalwateruse_h)
        # weekly profile
        for i, day in enumerate(
            [
                "monday",
                "tuesday",
                "wednesday",
                "thursday",
                "friday",
                "saturday",
                "sunday",
            ]
        ):
            set_df_value(f"daywatper", (i,), getattr(irrigation.daywatper, day))
            set_df_value(f"daywat", (i,), getattr(irrigation.daywat, day))
        # 24-hour profile
        for hour in range(24):
            for i, day in enumerate(["working_day", "holiday"]):
                for var in ["wuprofa_24hr", "wuprofm_24hr"]:
                    val = getattr(getattr(irrigation, var), day)[f"{hour+1}"]
                    set_df_value(var, (hour, i), val)

        # Process anthropogenic emissions
        anthro = props.anthropogenic_emissions
        for var in ["startdls", "enddls"]:
            set_df_value(var, 0, getattr(anthro, var))

        # heat
        anthro_heat = anthro.heat
        # day-specific mappings for anthro_heat
        daymap_anthro_heat = {
            "qf0_beu": anthro_heat.qf0_beu,
            "qf_a": anthro_heat.qf_a,
            "qf_b": anthro_heat.qf_b,
            "qf_c": anthro_heat.qf_c,
            "baset_cooling": anthro_heat.baset_cooling,
            "baset_heating": anthro_heat.baset_heating,
            "ah_min": anthro_heat.ah_min,
            "ah_slope_cooling": anthro_heat.ah_slope_cooling,
            "ah_slope_heating": anthro_heat.ah_slope_heating,
            "popdensdaytime": anthro_heat.popdensdaytime,
        }

        # Set values for working day and holiday
        for i, day in enumerate(["working_day", "holiday"]):
            for name, value in daymap_anthro_heat.items():
                set_df_value(name, (i,), getattr(value, day))

            # Hourly profiles
            for hour in range(24):
                set_df_value(
                    "ahprof_24hr",
                    (hour, i),
                    getattr(anthro_heat.ahprof_24hr, day)[f"{hour+1}"],
                )

        # CO2 parameters
        anthro_co2 = anthro.co2
        # invariant parameters
        map_anthro_co2 = {
            "co2pointsource": anthro_co2.co2pointsource,
            "ef_umolco2perj": anthro_co2.ef_umolco2perj,
            "enef_v_jkm": anthro_co2.enef_v_jkm,
            "trafficunits": anthro_co2.trafficunits,
            "frfossilfuel_heat": anthro_co2.frfossilfuel_heat,
            "frfossilfuel_nonheat": anthro_co2.frfossilfuel_nonheat,
            "maxfcmetab": anthro_co2.maxfcmetab,
            "maxqfmetab": anthro_co2.maxqfmetab,
            "minfcmetab": anthro_co2.minfcmetab,
            "minqfmetab": anthro_co2.minqfmetab,
        }
        for name, value in map_anthro_co2.items():
            set_df_value(name, 0, value)

        # Traffic and emission factors for weekday/weekend
        daymap_anthro_co2_traffic = {
            "trafficrate": anthro_co2.trafficrate,
            "fcef_v_kgkm": anthro_co2.fcef_v_kgkm,
        }
        for i, day in enumerate(["working_day", "holiday"]):
            for name, value in daymap_anthro_co2_traffic.items():
                set_df_value(name, (i,), getattr(value, day))

        # 24-hour profiles dependent on day type
        dayhourmap_anthro_co2_traffic = {
            "traffprof_24hr": anthro_co2.traffprof_24hr,
            "humactivity_24hr": anthro_co2.humactivity_24hr,
        }
        for i, day in enumerate(["working_day", "holiday"]):
            for hour in range(24):
                for name, value in dayhourmap_anthro_co2_traffic.items():
                    set_df_value(name, (hour, i), getattr(value, day)[f"{hour+1}"])

        # Snow parameters
        snow = props.snow
        map_snow = {
            "crwmax": snow.crwmax,
            "crwmin": snow.crwmin,
            "narp_emis_snow": snow.narp_emis_snow,
            "preciplimit": snow.preciplimit,
            "preciplimitalb": snow.preciplimitalb,
            "snowalbmax": snow.snowalbmax,
            "snowalbmin": snow.snowalbmin,
            "snowdensmin": snow.snowdensmin,
            "snowdensmax": snow.snowdensmax,
            "snowlimbldg": snow.snowlimbldg,
            "snowlimpaved": snow.snowlimpaved,
            "tau_a": snow.tau_a,
            "tau_r": snow.tau_r,
            "tau_f": snow.tau_f,
            "tempmeltfact": snow.tempmeltfact,
            "radmeltfact": snow.radmeltfact,
        }
        for name, value in map_snow.items():
            set_df_value(name, 0, value)

        # Missing height parameters for vegetation and buildings
        map_veg_height = {
            "bldgh": props.land_cover.bldgs.bldgh,  # Building height
            "evetreeh": props.land_cover.evetr.evetreeh,  # Evergreen tree height
            "dectreeh": props.land_cover.dectr.dectreeh,  # Deciduous tree height
        }
        for name, value in map_veg_height.items():
            set_df_value(name, 0, value)

        # Missing FAI parameters
        map_faibldg = {
            "faibldg": props.land_cover.bldgs.faibldg,
            "faievetree": props.land_cover.evetr.faievetree,
            "faidectree": props.land_cover.dectr.faidectree,
        }
        for name, value in map_faibldg.items():
            set_df_value(name, 0, value)

        # Population profile parameters
        set_df_value("popdensnighttime", 0, anthro_heat.popdensnighttime)
        for i, day in enumerate(["working_day", "holiday"]):
            # 24-hour population profile
            for hour in range(24):
                set_df_value(
                    "popprof_24hr",
                    (hour, i),
                    getattr(anthro_heat.popprof_24hr, day)[f"{hour+1}"],
                )

        dayhourmap_snow = {
            "snowprof_24hr": snow.snowprof_24hr,
        }
        for i, day in enumerate(["working_day", "holiday"]):
            for hour in range(24):
                for name, value in dayhourmap_snow.items():
                    set_df_value(name, (hour, i), getattr(value, day)[f"{hour+1}"])

        # Surface-specific parameter mappings
        def get_surface_mappings(surface, surf_idx):
            """Get mappings for surface parameters"""
            base_map = {
                "sfr_surf": ((surf_idx,), surface.sfr),
                "emis": ((surf_idx,), surface.emis),
            }

            # Optional parameters with direct mapping
            optional_params = {
                "alb": "alb",
                "ohm_threshsw": "ohm_threshsw",
                "ohm_threshwd": "ohm_threshwd",
                "chanohm": "chanohm",
                "cpanohm": "cpanohm",
                "kkanohm": "kkanohm",
                "soildepth": "soildepth",
                "soilstorecap": "soilstorecap_surf",
                "statelimit": "statelimit_surf",
                "wetthresh": "wetthresh_surf",
                "flowchange": "flowchange",
                "snowpacklimit": "snowpacklimit",
            }

            for attr, df_name in optional_params.items():
                if hasattr(surface, attr):
                    # special case for flowchange
                    if attr == "flowchange":
                        base_map[df_name] = (0, getattr(surface, attr))
                    else:
                        base_map[df_name] = ((surf_idx,), getattr(surface, attr))

            return base_map

        def get_vegetation_mappings(surface, surf_idx):
            """Get mappings specific to vegetation surfaces"""
            veg_map = {}

            # LAI parameters mapping
            idx = surf_idx - 2
            if hasattr(surface, "lai"):
                lai = surface.lai
                lai_params = {
                    "baset": lai.baset,
                    "gddfull": lai.gddfull,
                    "basete": lai.basete,
                    "sddfull": lai.sddfull,
                    "laimin": lai.laimin,
                    "laimax": lai.laimax,
                    "laitype": lai.laitype,
                }
                for param, value in lai_params.items():
                    veg_map[param] = ((idx,), value)

                # LAI power parameters
                lai_power_vars = [
                    "growth_lai",
                    "growth_gdd",
                    "senescence_lai",
                    "senescence_sdd",
                ]
                for i, var in enumerate(lai_power_vars):
                    veg_map["laipower"] = ((i, idx), getattr(lai.laipower, var))

            # CO2 parameters mapping
            co2_params = {
                "beta_bioco2": "beta_bioco2",
                "beta_enh_bioco2": "beta_enh_bioco2",
                "alpha_bioco2": "alpha_bioco2",
                "alpha_enh_bioco2": "alpha_enh_bioco2",
                "resp_a": "resp_a",
                "resp_b": "resp_b",
                "theta_bioco2": "theta_bioco2",
            }

            for attr, df_name in co2_params.items():
                if hasattr(surface, attr):
                    veg_map[df_name] = ((idx,), getattr(surface, attr))

            if hasattr(surface, "alb_min") and hasattr(surface, "alb_max"):
                # set_df_value("alb", (surf_idx,), surface.alb_min)  # Use min as default
                veg_map["alb"] = ((surf_idx,), surface.alb_min)
                surf_name = surface._surface_type.value
                str_alb_max = f"albmax_{surf_name}"
                str_alb_min = f"albmin_{surf_name}"
                veg_map[str_alb_max] = (0, surface.alb_max)
                veg_map[str_alb_min] = (0, surface.alb_min)

            return veg_map

        # Surface properties
        surface_map = {
            "paved": 0,
            "bldgs": 1,
            "evetr": 2,
            "dectr": 3,
            "grass": 4,
            "bsoil": 5,
            "water": 6,
        }

        # Process each surface type
        for surf_name, surf_idx in surface_map.items():
            surface = getattr(props.land_cover, surf_name)

            # # Get and apply base surface mappings
            surface_maps = get_surface_mappings(surface, surf_idx)

            for name, (indices, value) in surface_maps.items():
                set_df_value(name, indices, value)

            # a dummy row for extra surface ! not used but essential for the structural completeness
            set_df_value("ohm_threshsw", (7,), surface.ohm_threshsw)

            # For vegetation surfaces, apply additional mappings
            if surf_name in ["dectr", "evetr", "grass"]:
                veg_maps = get_vegetation_mappings(surface, surf_idx)
                for name, (indices, value) in veg_maps.items():
                    set_df_value(name, indices, value)

            # Add water distribution parameters
            if surface.waterdist is not None:
                # Map of possible distribution targets based on (8,6) dimensionality
                # 8 rows: Paved(0), Bldgs(1), EveTr(2), DecTr(3), Grass(4), BSoil(5), Water(6), Extra(7)
                # 6 cols: Paved(0), Bldgs(1), EveTr(2), DecTr(3), Grass(4), BSoil(5)
                dist_targets = {
                    "to_paved": 0,
                    "to_bldgs": 1,
                    "to_evetr": 2,
                    "to_dectr": 3,
                    "to_grass": 4,
                    "to_bsoil": 5,
                    "to_water": 6,
                }

                # Set water distribution values
                for target, value in surface.waterdist.__dict__.items():
                    try:
                        target_idx = dist_targets[target]
                        # print(surf_idx, target_idx, value)
                    except KeyError:
                        print(f"Target {target} not found in dist_targets")
                        continue
                    if value is not None:
                        if target == "to_runoff":
                            # Special case for runoff (separate column)
                            # set_df_value("waterdist_runoff", (surf_idx,), value)
                            pass
                        elif target == "to_soilstore":
                            # Special case for soil store (separate column)
                            # set_df_value("waterdist_soilstore", (surf_idx,), value)
                            pass
                        elif target == "to_water":
                            # Special case for water (row 6)
                            set_df_value("waterdist", (6, surf_idx), value)
                        elif target in dist_targets:
                            # Regular surface-to-surface distributions
                            target_idx = dist_targets[target]
                            # Note: in the DataFrame, first index is FROM surface, second is TO surface
                            set_df_value("waterdist", (target_idx, surf_idx), value)

                # Set unused row (7) to 0
                for col_idx in range(6):
                    set_df_value("waterdist", (7, col_idx), 0.0)
                # Set diagonal elements (targets to themselves) to 0
                for row_idx in range(6):
                    set_df_value("waterdist", (row_idx, row_idx), 0.0)

            # porosity parameters for deciduous trees
            if hasattr(surface, "pormin_dec"):
                set_df_value("pormin_dec", 0, surface.pormin_dec)
            if hasattr(surface, "pormax_dec"):
                set_df_value("pormax_dec", 0, surface.pormax_dec)

            # capacity parameters for deciduous trees
            if hasattr(surface, "capmin_dec"):
                set_df_value("capmin_dec", 0, surface.capmin_dec)
            if hasattr(surface, "capmax_dec"):
                set_df_value("capmax_dec", 0, surface.capmax_dec)
            if hasattr(surface, "min_res_bioco2"):
                set_df_value("min_res_bioco2", (surf_idx - 2,), surface.min_res_bioco2)

            # OHM coefficients
            if surface.ohm_coef:
                for i, (a1, a2, a3) in enumerate(
                    zip(
                        surface.ohm_coef.a1.values(),
                        surface.ohm_coef.a2.values(),
                        surface.ohm_coef.a3.values(),
                    )
                ):
                    set_df_value("ohm_coef", (surf_idx, i, 0), a1)
                    set_df_value("ohm_coef", (surf_idx, i, 1), a2)
                    set_df_value("ohm_coef", (surf_idx, i, 2), a3)
                    # dummy row for extra surface
                    set_df_value("ohm_coef", (7, i, 0), a1)
                    set_df_value("ohm_coef", (7, i, 1), a2)
                    set_df_value("ohm_coef", (7, i, 2), a3)

            # Storage and drain parameters
            if hasattr(surface, "storedrainprm"):
                for i, var in enumerate(
                    [
                        "store_min",
                        "store_max",
                        "store_cap",
                        "drain_eq",
                        "drain_coef_1",
                        "drain_coef_2",
                    ]
                ):
                    set_df_value(
                        "storedrainprm",
                        (i, surf_idx),
                        getattr(surface.storedrainprm, var),
                    )
            # Irrigation coefficients
            if hasattr(surface, "ie_a"):
                set_df_value("ie_a", (surf_idx - 2,), surface.ie_a)
            if hasattr(surface, "ie_m"):
                set_df_value("ie_m", (surf_idx - 2,), surface.ie_m)

            # Maximum conductance
            if hasattr(surface, "maxconductance"):
                idx = surf_idx - 2
                set_df_value("maxconductance", (idx,), surface.maxconductance)

            # Snow pack limit
            if hasattr(surface, "snowpacklimit"):
                set_df_value("snowpacklimit", (surf_idx,), surface.snowpacklimit)

            # LAI parameters
            if hasattr(surface, "lai"):
                lai = surface.lai
                idx = surf_idx - 2
                set_df_value("baset", (idx,), lai.baset)
                set_df_value("gddfull", (idx,), lai.gddfull)
                set_df_value("basete", (idx,), lai.basete)
                set_df_value("sddfull", (idx,), lai.sddfull)
                set_df_value("laimin", (idx,), lai.laimin)
                set_df_value("laimax", (idx,), lai.laimax)
                set_df_value("laitype", (idx,), lai.laitype)
                for i, var in enumerate(
                    [
                        "growth_lai",
                        "growth_gdd",
                        "senescence_lai",
                        "senescence_sdd",
                    ]
                ):
                    set_df_value("laipower", (i, idx), getattr(lai.laipower, var))

            # CO2 parameters for vegetated surfaces
            if hasattr(surface, "beta_bioco2"):
                idx = surf_idx - 2
                set_df_value("beta_bioco2", (idx,), surface.beta_bioco2)
                set_df_value("beta_enh_bioco2", (idx,), surface.beta_enh_bioco2)
                set_df_value("alpha_bioco2", (idx,), surface.alpha_bioco2)
                set_df_value("alpha_enh_bioco2", (idx,), surface.alpha_enh_bioco2)

            # Add to surface properties section
            if hasattr(surface, "sathydraulicconduct"):
                set_df_value(
                    "sathydraulicconduct", (surf_idx,), surface.sathydraulicconduct
                )

            # Water specific parameters
            if hasattr(surface, "flowchange"):
                set_df_value("flowchange", 0, surface.flowchange)

            # Add to CO2 parameters section for vegetated surfaces
            if hasattr(surface, "resp_a"):
                set_df_value("resp_a", (surf_idx - 2,), surface.resp_a)
            if hasattr(surface, "resp_b"):
                set_df_value("resp_b", (surf_idx - 2,), surface.resp_b)
            if hasattr(surface, "theta_bioco2"):
                set_df_value("theta_bioco2", (surf_idx - 2,), surface.theta_bioco2)

            # Snowpack limit
            if hasattr(surface, "snowpacklimit"):
                set_df_value("snowpacklimit", (surf_idx,), surface.snowpacklimit)

            # irrigation fraction
            if hasattr(surface, "irrfrac"):
                set_df_value(f"irrfrac{surf_name}", 0, surface.irrfrac)

            # Thermal layers
            if hasattr(surface, "thermal_layers"):
                # Get the lists from thermal_layers
                dz_list = surface.thermal_layers.dz
                k_list = surface.thermal_layers.k
                cp_list = surface.thermal_layers.cp

                # Set each value in the lists
                for i in range(5):  # We know there are exactly 5 values
                    set_df_value("dz_surf", (surf_idx, i), dz_list[i])
                    set_df_value("k_surf", (surf_idx, i), k_list[i])
                    set_df_value("cp_surf", (surf_idx, i), cp_list[i])

            # Initial states
            init_state = getattr(site.initial_states, surf_name)
            if init_state:
                # Basic state parameters
                set_df_value("state_surf", (surf_idx,), init_state.state)
                set_df_value("soilstore_surf", (surf_idx,), init_state.soilstore)

                # fill in dummy values for variables that are not needed for users
                for i in range(12):
                    set_df_value("hdd_id", (i,), 0)
                set_df_value("dqndt", 0, 1)
                set_df_value("dqnsdt", 0, 10)
                set_df_value("qn_av", 0, 10)
                set_df_value("qn_s_av", 0, 10)
                set_df_value("dt_since_start", 0, 0)
                set_df_value("lenday_id", 0, 0)
                set_df_value("tmax_id", 0, 0)
                set_df_value("tmin_id", 0, 0)
                set_df_value("tstep_prev", 0, 300)
                set_df_value("tair_av", 0, 0)
                set_df_value("snowfallcum", 0, 0)

                # Snow-related parameters
                if init_state.snowfrac is not None:
                    set_df_value("snowfrac", (surf_idx,), init_state.snowfrac)
                if init_state.snowpack is not None:
                    set_df_value("snowpack", (surf_idx,), init_state.snowpack)
                if init_state.snowwater is not None:
                    set_df_value("snowwater", (surf_idx,), init_state.snowwater)
                if init_state.snowdens is not None:
                    set_df_value("snowdens", (surf_idx,), init_state.snowdens)
                if init_state.icefrac is not None:
                    set_df_value("icefrac", (surf_idx,), init_state.icefrac)

                # Vegetation-specific parameters
                if surf_name in ["dectr", "evetr", "grass"]:
                    # albedo
                    if init_state.alb_id is not None:
                        set_df_value(f"alb{surf_name}_id", 0, init_state.alb_id)
                    # LAI
                    if init_state.lai_id is not None:
                        set_df_value(f"lai_id", (surf_idx - 2,), init_state.lai_id)
                    # GDD
                    if init_state.gdd_id is not None:
                        set_df_value(f"gdd_id", (surf_idx - 2,), init_state.gdd_id)
                    # SDD
                    if init_state.sdd_id is not None:
                        set_df_value(f"sdd_id", (surf_idx - 2,), init_state.sdd_id)
                    # water use
                    if init_state.wu is not None:
                        set_df_value(
                            f"wuday_id",
                            ((surf_idx - 2) * 3 + 0,),
                            init_state.wu.wu_total,
                        )
                        set_df_value(
                            f"wuday_id",
                            ((surf_idx - 2) * 3 + 1,),
                            init_state.wu.wu_auto,
                        )
                        set_df_value(
                            f"wuday_id",
                            ((surf_idx - 2) * 3 + 2,),
                            init_state.wu.wu_manual,
                        )
                    # Additional parameters for deciduous trees
                    if surf_name == "dectr":
                        if init_state.decidcap_id is not None:
                            set_df_value("decidcap_id", 0, init_state.decidcap_id)
                        if init_state.porosity_id is not None:
                            set_df_value("porosity_id", 0, init_state.porosity_id)
                # temperature
                for k, var in enumerate(["temperature", "tsfc", "tin"]):
                    if var == "tsfc":
                        set_df_value(f"tsfc_surf", (surf_idx,), init_state.tsfc)
                    elif var == "temperature":
                        for k in range(5):
                            set_df_value(
                                f"temp_surf", (surf_idx, k), getattr(init_state, var)[k]
                            )
                    elif var == "tin":
                        set_df_value(f"tin_surf", (surf_idx,), init_state.tin)
            # Set initial snow albedo
            set_df_value("snowalb", 0, site.initial_states.snowalb)

        # Vertical layers
        if hasattr(props, "vertical_layers"):
            print("vertical_layers here")
            vertical_layers = props.vertical_layers
            set_df_value("nlayer", 0, vertical_layers.nlayer)
            for i in range(vertical_layers.nlayer + 1):
                set_df_value("height", (i,), vertical_layers.height[i])
            for i in range(vertical_layers.nlayer):
                set_df_value("building_scale", (i,), vertical_layers.building_scale[i])
                set_df_value("building_frac", (i,), vertical_layers.building_frac[i])
                set_df_value("veg_scale", (i,), vertical_layers.veg_scale[i])
                set_df_value("veg_frac", (i,), vertical_layers.veg_frac[i])
            for i, layer in enumerate(vertical_layers.roofs):
                set_df_value(
                    f"roof_albedo_dir_mult_fact",
                    (0, i),
                    layer.roof_albedo_dir_mult_fact,
                )
                set_df_value(f"alb_roof", (i,), layer.alb)
                set_df_value(f"emis_roof", (i,), layer.emis)
                set_df_value(f"statelimit_roof", (i,), layer.statelimit)
                set_df_value(f"soilstorecap_roof", (i,), layer.soilstorecap)
                set_df_value(f"wetthresh_roof", (i,), layer.wetthresh)
                thermal_layers = layer.thermal_layers
                for j, var in enumerate(["dz", "k", "cp"]):
                    for k in range(5):
                        set_df_value(
                            f"{var}_roof", (i, k), getattr(thermal_layers, var)[k]
                        )
            for i, layer in enumerate(vertical_layers.walls):
                set_df_value(f"wall_specular_frac", (0, i), layer.wall_specular_frac)
                set_df_value(f"alb_wall", (i,), layer.alb)
                set_df_value(f"emis_wall", (i,), layer.emis)
                set_df_value(f"statelimit_wall", (i,), layer.statelimit)
                set_df_value(f"soilstorecap_wall", (i,), layer.soilstorecap)
                set_df_value(f"wetthresh_wall", (i,), layer.wetthresh)
                thermal_layers = layer.thermal_layers
                for j, var in enumerate(["dz", "k", "cp"]):
                    for k in range(5):
                        set_df_value(
                            f"{var}_wall", (i, k), getattr(thermal_layers, var)[k]
                        )
        for building_facet in ["roofs", "walls"]:
            facet = building_facet[:-1]
            if hasattr(site.initial_states, building_facet):
                for j, layer in enumerate(getattr(site.initial_states, building_facet)):
                    set_df_value(f"state_{facet}", (j,), layer.state)
                    set_df_value(f"soilstore_{facet}", (j,), layer.soilstore)
                    # set_df_value(f"snowwater_{facet}", (j,), layer.snowwater)
                    # set_df_value(f"snowdens_{facet}", (j,), layer.snowdens)
                    # set_df_value(f"snowfrac_{facet}", (j,), layer.snowfrac)
                    # set_df_value(f"snowpack_{facet}", (j,), layer.snowpack)
                    for k, var in enumerate(["temperature", "tsfc", "tin"]):
                        if var == "tsfc":
                            set_df_value(f"tsfc_{facet}", (j,), layer.tsfc)
                        elif var == "temperature":
                            for k in range(5):
                                set_df_value(
                                    f"temp_{facet}", (j, k), getattr(layer, var)[k]
                                )
                        elif var == "tin":
                            set_df_value(f"tin_{facet}", (j,), layer.tin)
        return df

    def to_df_state_new(self) -> pd.DataFrame:
        """Convert config to DataFrame state format"""
        list_df_site = []
        for grid_id in range(len(self.site)):
            df_site = self.site[grid_id].to_df_state(grid_id)
            list_df_site.append(df_site)
        df = pd.concat(list_df_site, axis=1)
        return df

    @classmethod
    def from_df_state(cls, df: pd.DataFrame) -> "SUEWSConfig":
        """Create config from DataFrame state"""
        # TODO: add from_df_state
        pass


if __name__ == "__main__":
    # Create list for collecting all exceptions
    exceptions = []

    # test the sample config
    # Load YAML config
    with open("./config-suews.yml", "r") as file:
        yaml_config = yaml.safe_load(file)

    # Create SUEWSConfig object
    suews_config = SUEWSConfig(**yaml_config[0])

    if exceptions:
        raise ExceptionGroup("Validation errors occurred", exceptions)

    print(r"testing suews_config done!")

    # pdb.set_trace()

    # Convert to DataFrame
    df_state_test = suews_config.to_df_state()
    df_state_test.to_pickle("./df_state_test.pkl")
    print("testing df_state done!")

    # checking if all properties are properly converted
    # Get the column differences
    df_state = pd.read_pickle("./df_state.pkl")
    df_state_cols = set(df_state.columns)
    df_test_cols = set(df_state_test.columns)

    print("Columns only in df_state:")
    print(sorted(df_state_cols - df_test_cols))

    print("\nColumns only in df_state_test:")
    print(sorted(df_test_cols - df_state_cols))

    print("\nTotal columns in df_state:", len(df_state_cols))
    print("Total columns in df_state_test:", len(df_test_cols))

    # Get columns with any NA values
    na_cols = df_state_test.columns[df_state_test.isna().any()].tolist()

    # Sort and print the column names
    print(f"{len(na_cols)} Columns containing NA values:")
    for col in sorted(na_cols):
        print(col)

    # test running supy
    import supy as sp

    df_state, df_forcing = sp.load_SampleData()
    df_state_test = pd.read_pickle("./df_state_test.pkl")
    sp.run_supy(df_forcing.iloc[0 : 288 * 10], df_state_test)


# # Convert back to config
# suews_config_back = SUEWSConfig.from_df_state(df_state)
# print("testing from_df_state done!")
