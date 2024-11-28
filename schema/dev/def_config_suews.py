from typing import Dict, List, Optional, Union, Literal, Tuple, Type
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
    col = pd.MultiIndex.from_tuples([("gridiv", "0")], names=["var", "ind_dim"])
    df_state = pd.DataFrame(index=idx, columns=col)
    df_state.loc[grid_id, ("gridiv", "0")] = grid_id
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
        return cls(snowalb=snowalb)


class WaterUse(BaseModel):
    wu_total: float = Field(ge=0, description="Total water use", default=0.0)
    wu_auto: float = Field(ge=0, description="Automatic water use", default=0.0)
    wu_manual: float = Field(ge=0, description="Manual water use", default=0.0)

    def to_df_state(self, veg_idx: int, grid_id: int) -> pd.DataFrame:
        """Convert water use to DataFrame state format."""
        df_state = init_df_state(grid_id)
        df_state.loc[grid_id, ("wuday_id", f"({veg_idx * 3 + 0},)")] = self.wu_total
        df_state.loc[grid_id, ("wuday_id", f"({veg_idx * 3 + 1},)")] = self.wu_auto
        df_state.loc[grid_id, ("wuday_id", f"({veg_idx * 3 + 2},)")] = self.wu_manual
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
        wu_total = df.loc[grid_id, ("wuday_id", f"({veg_idx * 3 + 0},)")]
        wu_auto = df.loc[grid_id, ("wuday_id", f"({veg_idx * 3 + 1},)")]
        wu_manual = df.loc[grid_id, ("wuday_id", f"({veg_idx * 3 + 2},)")]

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

    def to_df_state(
        self, grid_id: int, vert_idx: int = None, is_roof: bool = False
    ) -> pd.DataFrame:
        """Convert base surface initial state to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing initial state parameters
        """
        df_state = init_df_state(grid_id)

        # Get surface index
        if vert_idx is None:
            idx = self.get_surface_index()
            str_type = "surf"
        else:
            idx = vert_idx
            str_type = "roof" if is_roof else "wall"
        # Set basic state parameters
        df_state[(f"state_{str_type}", f"({idx},)")] = self.state
        df_state[(f"soilstore_{str_type}", f"({idx},)")] = self.soilstore

        # Set snow/ice parameters if present
        if self.snowfrac is not None:
            df_state[(f"snowfrac", f"({idx},)")] = self.snowfrac
        if self.snowpack is not None:
            df_state[(f"snowpack", f"({idx},)")] = self.snowpack
        if self.icefrac is not None:
            df_state[(f"icefrac", f"({idx},)")] = self.icefrac
        if self.snowwater is not None:
            df_state[(f"snowwater", f"({idx},)")] = self.snowwater
        if self.snowdens is not None:
            df_state[(f"snowdens", f"({idx},)")] = self.snowdens

        # Set temperature parameters
        for i, temp in enumerate(self.temperature):
            df_state[(f"temp_{str_type}", f"({idx}, {i})")] = temp

        if self.tsfc is not None:
            df_state[(f"tsfc_{str_type}", f"({idx},)")] = self.tsfc
        if self.tin is not None:
            df_state[(f"tin_{str_type}", f"({idx},)")] = self.tin

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int) -> "SurfaceInitialState":
        """
        Reconstruct SurfaceInitialState from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing surface state parameters.
            grid_id (int): Grid ID for the DataFrame index.
            surf_idx (int): Surface index for identifying columns.
            str_type (str): Surface type prefix ("surf", "roof", or "wall").

        Returns:
            SurfaceInitialState: Instance of SurfaceInitialState.
        """
        # Base surface state parameters
        state = df.loc[grid_id, ("state_surf", f"({surf_idx},)")]
        soilstore = df.loc[grid_id, ("soilstore_surf", f"({surf_idx},)")]

        # Snow/ice parameters
        snowfrac = df.loc[grid_id, ("snowfrac", f"({surf_idx},)")]
        snowpack = df.loc[grid_id, ("snowpack", f"({surf_idx},)")]
        icefrac = df.loc[grid_id, ("icefrac", f"({surf_idx},)")]
        snowwater = df.loc[grid_id, ("snowwater", f"({surf_idx},)")]
        snowdens = df.loc[grid_id, ("snowdens", f"({surf_idx},)")]

        # Temperature parameters
        temperature = [
            df.loc[grid_id, ("temp_surf", f"({surf_idx}, {i})")] for i in range(5)
        ]

        # Exterior and interior surface temperature
        tsfc = df.loc[grid_id, ("tsfc_surf", f"({surf_idx},)")]
        tin = df.loc[grid_id, ("tin_surf", f"({surf_idx},)")]

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

class InitialStatePaved(SurfaceInitialState):
    _surface_type: Literal[SurfaceType.PAVED] = SurfaceType.PAVED


class InitialStateBldgs(SurfaceInitialState):
    _surface_type: Literal[SurfaceType.BLDGS] = SurfaceType.BLDGS


class VegInitialState(SurfaceInitialState):
    """Base initial state parameters for vegetated surfaces"""

    alb_id: float = Field(
        description="Initial albedo for vegetated surfaces", default=0.1
    )
    lai_id: float = Field(description="Initial leaf area index", default=1.0)
    gdd_id: float = Field(description="Growing degree days ID", default=0)
    sdd_id: float = Field(description="Senescence degree days ID", default=0)
    wu: WaterUse = Field(default_factory=WaterUse)

    @model_validator(mode="after")
    def validate_surface_state(self) -> "VegInitialState":
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
        veg_idx = surf_idx - 2

        # Add vegetated surface specific parameters
        # alb is universal so use surf_idx
        df_state[("alb", f"({surf_idx},)")] = self.alb_id
        # others are aligned with veg_idx
        df_state[("lai_id", f"({veg_idx},)")] = self.lai_id
        df_state[("gdd_id", f"({veg_idx},)")] = self.gdd_id
        df_state[("sdd_id", f"({veg_idx},)")] = self.sdd_id

        # Add water use parameters
        df_wu = self.wu.to_df_state(veg_idx, grid_id)
        df_state = pd.concat([df_state, df_wu], axis=1)

        # Drop any duplicate columns
        df_state = df_state.loc[:, ~df_state.columns.duplicated()]

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int
    ) -> "VegInitialState":
        """
        Reconstruct VegetatedSurfaceInitialState from a DataFrame state format."""
        base_instance = SurfaceInitialState.from_df_state(df, grid_id, surf_idx)

        # Vegetated surface-specific parameters
        alb_key = ("alb", f"({surf_idx},)")
        lai_key = ("lai_id", f"({surf_idx - 2},)")
        gdd_key = ("gdd_id", f"({surf_idx - 2},)")
        sdd_key = ("sdd_id", f"({surf_idx - 2},)")

        alb_id = df.loc[grid_id, alb_key]
        lai_id = df.loc[grid_id, lai_key]
        gdd_id = df.loc[grid_id, gdd_key]
        sdd_id = df.loc[grid_id, sdd_key]

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


class InitialStateEvetr(VegInitialState):
    _surface_type: Literal[SurfaceType.EVETR] = SurfaceType.EVETR

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert evergreen tree initial state to DataFrame state format."""
        df_state = super().to_df_state(grid_id)
        df_state[("albevetr_id", "0")] = self.alb_id
        return df_state
    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int, surf_idx: int) -> "InitialStateEvetr":
        """
        Reconstruct InitialStateEvetr from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing EVETR state parameters.
            grid_id (int): Grid ID for the DataFrame index.
            surf_idx (int): Surface index for identifying columns.

        Returns:
            InitialStateEvetr: Instance of InitialStateEvetr.
        """
        # Call the parent class to extract common fields
        base_instance = super().from_df_state(df, grid_id, surf_idx)

        # Extract the EVETR-specific field
        alb_id = df.loc[grid_id, ("albevetr_id", "0")].item()

        # Use `base_instance.dict()` to pass the existing attributes, excluding `alb_id` to avoid duplication
        base_instance_dict = base_instance.dict()
        base_instance_dict["alb_id"] = alb_id  # Update alb_id explicitly

        # Return a new instance with the updated dictionary
        return cls(**base_instance_dict)


class InitialStateDectr(VegInitialState):
    """Initial state parameters for deciduous trees"""

    porosity_id: float = Field(
        default=0.2, description="Initial porosity for deciduous trees"
    )
    decidcap_id: float = Field(
        default=0.3, description="Initial deciduous capacity for deciduous trees"
    )
    _surface_type: Literal[SurfaceType.DECTR] = SurfaceType.DECTR

    @model_validator(mode="after")
    def validate_surface_state(self) -> "InitialStateDectr":
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

        # Add deciduous tree specific parameters
        df_state[("porosity_id", "0")] = self.porosity_id
        df_state[("decidcap_id", "0")] = self.decidcap_id
        df_state[("albdectr_id", "0")] = self.alb_id

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int
    ) -> "InitialStateDectr":
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
        base_instance = VegInitialState.from_df_state(df, grid_id, surf_idx)

        # Deciduous tree-specific parameters
        porosity_id = df.loc[grid_id, ("porosity_id", "0")]
        decidcap_id = df.loc[grid_id, ("decidcap_id", "0")]

        return cls(
            **base_instance.dict(),
            porosity_id=porosity_id,
            decidcap_id=decidcap_id,
        )


class InitialStateGrass(VegInitialState):
    _surface_type: Literal[SurfaceType.GRASS] = SurfaceType.GRASS

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert grass initial state to DataFrame state format."""
        df_state = super().to_df_state(grid_id)
        df_state[("albgrass_id", "0")] = self.alb_id
        return df_state
    
    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int, surf_idx: int) -> "InitialStateGrass":
        """
        Reconstruct InitialStateGrass from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing grass state parameters.
            grid_id (int): Grid ID for the DataFrame index.
            surf_idx (int): Surface index for identifying columns.

        Returns:
            InitialStateGrass: Instance of InitialStateGrass.
        """
        # Call the parent class to extract common fields
        base_instance = super().from_df_state(df, grid_id, surf_idx)

        # Extract the GRASS-specific field
        alb_id = df.loc[grid_id, ("albgrass_id", "0")].item()

        # Use `base_instance.dict()` to pass the existing attributes, excluding `alb_id` to avoid duplication
        base_instance_dict = base_instance.dict()
        base_instance_dict["alb_id"] = alb_id  # Update alb_id explicitly

        # Return a new instance with the updated dictionary
        return cls(**base_instance_dict)



class InitialStateBsoil(SurfaceInitialState):
    _surface_type: Literal[SurfaceType.BSOIL] = SurfaceType.BSOIL


class InitialStateWater(SurfaceInitialState):
    _surface_type: Literal[SurfaceType.WATER] = SurfaceType.WATER


class InitialStates(BaseModel):
    """Initial conditions for the SUEWS model"""

    snowalb: float = Field(ge=0, le=1, description="Initial snow albedo", default=0.5)
    paved: InitialStatePaved = Field(default_factory=InitialStatePaved)
    bldgs: InitialStateBldgs = Field(default_factory=InitialStateBldgs)
    evetr: InitialStateEvetr = Field(default_factory=InitialStateEvetr)
    dectr: InitialStateDectr = Field(default_factory=InitialStateDectr)
    grass: InitialStateGrass = Field(default_factory=InitialStateGrass)
    bsoil: InitialStateBsoil = Field(default_factory=InitialStateBsoil)
    water: InitialStateWater = Field(default_factory=InitialStateWater)
    roofs: Optional[List[SurfaceInitialState]] = Field(
        default=[
            SurfaceInitialState(),
            SurfaceInitialState(),
            SurfaceInitialState(),
        ],
        description="Initial states for roof layers",
    )
    walls: Optional[List[SurfaceInitialState]] = Field(
        default=[
            SurfaceInitialState(),
            SurfaceInitialState(),
            SurfaceInitialState(),
        ],
        description="Initial states for wall layers",
    )

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
        for facet_list, facet_type in [(self.roofs, "roof"), (self.walls, "wall")]:
            if facet_list is not None:  # Check for None explicitly
                for i, facet in enumerate(facet_list):
                    is_roof = facet_type == "roof"
                    df_facet = facet.to_df_state(grid_id, i, is_roof)
                    df_state = pd.concat([df_state, df_facet], axis=1)
                    df_state = df_state.sort_index(axis=1)

        # add dummy columns to conform to SUEWS convention
        list_cols = [
            "dqndt",
            "dqnsdt",
            "dt_since_start",
            "lenday_id",
            "qn_av",
            "qn_s_av",
            "tair_av",
            "tmax_id",
            "tmin_id",
            "tstep_prev",
            "snowfallcum",
        ]
        for col in list_cols:
            df_state[(col, "0")] = 0
            df_state = df_state.sort_index(axis=1)
        # special treatment for hdd_id
        for i in range(12):
            df_state[(f"hdd_id", f"({i},)")] = 0
            df_state = df_state.sort_index(axis=1)
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

        # Correct surface types to their respective classes
        surface_types = {
            "paved": InitialStatePaved,
            "bldgs": InitialStateBldgs,
            "evetr": InitialStateEvetr,
            "dectr": InitialStateDectr,
            "grass": InitialStateGrass,
            "bsoil": InitialStateBsoil,
            "water": InitialStateWater,
        }

        # Reconstruct surface states
        surfaces = {
            name: surface_class.from_df_state(df, grid_id, idx)
            for idx, (name, surface_class) in enumerate(surface_types.items())
        }

        # Reconstruct roof and wall states
        def reconstruct_layers(layer_name: str, surface_class: Type[SurfaceInitialState]):
            layers = []
            idx = 0
            while True:
                try:
                    layer = surface_class.from_df_state(df, grid_id, idx)
                    layers.append(layer)
                    idx += 1
                except KeyError:
                    break
            return layers

        roofs = reconstruct_layers("roof", SurfaceInitialState)
        walls = reconstruct_layers("wall", SurfaceInitialState)

        # Ensure each surface is passed to the corresponding field
        return cls(
            snowalb=snowalb,
            paved=surfaces["paved"],
            bldgs=surfaces["bldgs"],
            evetr=surfaces["evetr"],
            dectr=surfaces["dectr"],
            grass=surfaces["grass"],
            bsoil=surfaces["bsoil"],
            water=surfaces["water"],
            roofs=roofs,
            walls=walls,
        )


class ThermalLayers(BaseModel):
    dz: List[float] = Field([0.1, 0.2, 0.3, 0.4, 0.5], min_items=5, max_items=5)
    k: List[float] = Field([1.0, 1.0, 1.0, 1.0, 1.0], min_items=5, max_items=5)
    cp: List[float] = Field([1000, 1000, 1000, 1000, 1000], min_items=5, max_items=5)

    def to_df_state(
        self,
        grid_id: int,
        idx: int,
        surf_type: Literal[
            "paved",
            "bldgs",
            "evetr",
            "dectr",
            "grass",
            "bsoil",
            "water",
            "roof",
            "wall",
        ],
    ) -> pd.DataFrame:
        """Convert thermal layer parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            surf_type: Surface type or facet type ("roof" or "wall")

        Returns:
            pd.DataFrame: DataFrame containing thermal layer parameters
        """
        df_state = init_df_state(grid_id)

        if surf_type == "roof":
            suffix = "roof"
        elif surf_type == "wall":
            suffix = "wall"
        else:
            suffix = "surf"

        # Add thermal layer parameters
        for i in range(5):
            df_state[(f"dz_{suffix}", f"({idx}, {i})")] = self.dz[i]
            df_state[(f"k_{suffix}", f"({idx}, {i})")] = self.k[i]
            df_state[(f"cp_{suffix}", f"({idx}, {i})")] = self.cp[i]

        return df_state

    @classmethod
    def from_df_state(
        cls,
        df: pd.DataFrame,
        grid_id: int,
        idx: int,
        surf_type: Union[SurfaceType, Literal["roof", "wall"]],
    ) -> "ThermalLayers":
        """Reconstruct ThermalLayers instance from DataFrame.

        Args:
            df: DataFrame containing thermal layer parameters.
            grid_id: Grid ID for the DataFrame index.
            idx: Surface index for identifying columns.
            surf_type: Surface type or facet type ("roof" or "wall").

        Returns:
            ThermalLayers: Reconstructed ThermalLayers instance.
        """
        dz = []
        k = []
        cp = []

        # Determine suffix based on surf_type
        if surf_type == "roof":
            suffix = "roof"
        elif surf_type == "wall":
            suffix = "wall"
        else:
            suffix = "surf"

        # Extract thermal layer parameters for each of the 5 layers
        for i in range(5):
            dz.append(df.loc[grid_id, (f"dz_{suffix}", f"({idx}, {i})")])
            k.append(df.loc[grid_id, (f"k_{suffix}", f"({idx}, {i})")])
            cp.append(df.loc[grid_id, (f"cp_{suffix}", f"({idx}, {i})")])

        # Return reconstructed instance
        return cls(dz=dz, k=k, cp=cp)


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

    def __init__(self, surface_type: SurfaceType):
        # Default distributions based on surface type
        default_distributions = {
            SurfaceType.PAVED: {
                "to_bldgs": 0.2,
                "to_evetr": 0.1,
                "to_dectr": 0.1,
                "to_grass": 0.1,
                "to_bsoil": 0.1,
                "to_water": 0.1,
                "to_runoff": 0.3,
            },
            SurfaceType.BLDGS: {
                "to_paved": 0.2,
                "to_evetr": 0.1,
                "to_dectr": 0.1,
                "to_grass": 0.1,
                "to_bsoil": 0.1,
                "to_water": 0.1,
                "to_runoff": 0.3,
            },
            SurfaceType.EVETR: {
                "to_paved": 0.1,
                "to_bldgs": 0.1,
                "to_dectr": 0.1,
                "to_grass": 0.1,
                "to_bsoil": 0.1,
                "to_water": 0.1,
                "to_soilstore": 0.4,
            },
            SurfaceType.DECTR: {
                "to_paved": 0.1,
                "to_bldgs": 0.1,
                "to_evetr": 0.1,
                "to_grass": 0.1,
                "to_bsoil": 0.1,
                "to_water": 0.1,
                "to_soilstore": 0.4,
            },
            SurfaceType.GRASS: {
                "to_paved": 0.1,
                "to_bldgs": 0.1,
                "to_dectr": 0.1,
                "to_evetr": 0.1,
                "to_bsoil": 0.1,
                "to_water": 0.1,
                "to_soilstore": 0.4,
            },
            SurfaceType.BSOIL: {
                "to_paved": 0.1,
                "to_bldgs": 0.1,
                "to_dectr": 0.1,
                "to_evetr": 0.1,
                "to_grass": 0.1,
                "to_water": 0.1,
                "to_soilstore": 0.4,
            },
        }

        # If surface type is provided, use its default distribution
        # surface_type = data.get('_surface_type')
        if surface_type and surface_type in default_distributions:
            # Merge provided data with defaults, prioritising provided data
            merged_data = {**default_distributions[surface_type]}
            super().__init__(**merged_data)
        else:
            # If no surface type or invalid surface type, just use provided data
            super().__init__()

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
                # "to_soilstore",
                # "to_runoff",
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
                param_tuples.append(("waterdist", f"({i}, {surf_idx})"))
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
            ("storedrainprm", f"({i}, {surf_idx})")
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
            df.loc[grid_id, ("storedrainprm", f"({i}, {surf_idx})")] = getattr(
                self, var
            )

        return df

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int
    ) -> "StorageDrainParams":
        """
        Reconstruct StorageDrainParams from DataFrame state format.

        Args:
            df: DataFrame containing storage and drain parameters.
            grid_id: Grid ID for the DataFrame index.
            surf_idx: Surface index (0=paved, 1=bldgs, 2=dectr, 3=evetr, 4=grass, 5=bsoil, 6=water).

        Returns:
            StorageDrainParams: Instance of StorageDrainParams.
        """
        # Define the parameter names and their indices
        param_map = {
            "store_min": 0,
            "store_max": 1,
            "store_cap": 2,
            "drain_eq": 3,
            "drain_coef_1": 4,
            "drain_coef_2": 5,
        }

        # Extract the values from the DataFrame
        params = {
            param: df.loc[grid_id, ("storedrainprm", f"({idx}, {surf_idx})")]
            for param, idx in param_map.items()
        }

        # Create an instance using the extracted parameters
        return cls(**params)


class OHM_Coefficient_season_wetness(BaseModel):
    summer_dry: float = Field(
        default=0.0, description="OHM coefficient for summer dry conditions"
    )
    summer_wet: float = Field(
        default=0.0, description="OHM coefficient for summer wet conditions"
    )
    winter_dry: float = Field(
        default=0.0, description="OHM coefficient for winter dry conditions"
    )
    winter_wet: float = Field(
        default=0.0, description="OHM coefficient for winter wet conditions"
    )

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
            str_idx = f"({surf_idx}, {idx}, {idx_a})"
            df_state.loc[grid_id, ("ohm_coef", str_idx)] = getattr(self, season_wetdry)

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int, idx_a: int
    ) -> "OHM_Coefficient_season_wetness":
        """
        Reconstruct OHM_Coefficient_season_wetness from DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing OHM coefficients.
            grid_id (int): Grid ID.
            surf_idx (int): Surface index.
            idx_a (int): Index for coefficient (0=a1, 1=a2, 2=a3).

        Returns:
            OHM_Coefficient_season_wetness: Reconstructed instance.
        """
        season_wetdry_map = {
            "summer_dry": 0,
            "summer_wet": 1,
            "winter_dry": 2,
            "winter_wet": 3,
        }

        # Extract values for each season/wetness combination
        params = {
            season_wetdry: df.loc[
                grid_id, ("ohm_coef", f"({surf_idx}, {idx}, {idx_a})")
            ]
            for season_wetdry, idx in season_wetdry_map.items()
        }

        return cls(**params)


class OHMCoefficients(BaseModel):
    a1: OHM_Coefficient_season_wetness = Field(
        default_factory=OHM_Coefficient_season_wetness,
        description="OHM coefficient a1 for different seasons and wetness conditions",
    )
    a2: OHM_Coefficient_season_wetness = Field(
        default_factory=OHM_Coefficient_season_wetness,
        description="OHM coefficient a2 for different seasons and wetness conditions",
    )
    a3: OHM_Coefficient_season_wetness = Field(
        default_factory=OHM_Coefficient_season_wetness,
        description="OHM coefficient a3 for different seasons and wetness conditions",
    )

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
            df_coef_extra = coef.to_df_state(
                grid_id, 7, idx_a
            )  # always include this extra row to conform to SUEWS convention
            df_state = pd.concat([df_state, df_coef, df_coef_extra], axis=1)

        # drop duplicate columns
        df_state = df_state.loc[:, ~df_state.columns.duplicated()]

        return df_state

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int
    ) -> "OHMCoefficients":
        """
        Reconstruct OHMCoefficients from DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing OHM coefficients.
            grid_id (int): Grid ID.
            surf_idx (int): Surface index.

        Returns:
            OHMCoefficients: Reconstructed instance.
        """
        # Reconstruct each coefficient (a1, a2, a3)
        a1 = OHM_Coefficient_season_wetness.from_df_state(df, grid_id, surf_idx, 0)
        a2 = OHM_Coefficient_season_wetness.from_df_state(df, grid_id, surf_idx, 1)
        a3 = OHM_Coefficient_season_wetness.from_df_state(df, grid_id, surf_idx, 2)

        return cls(a1=a1, a2=a2, a3=a3)


class SurfaceProperties(BaseModel):
    """Base properties for all surface types"""

    sfr: float = Field(ge=0, le=1, description="Surface fraction", default=1.0 / 7)
    emis: float = Field(ge=0, le=1, description="Surface emissivity", default=0.95)
    chanohm: Optional[float] = Field(default=0.0)
    cpanohm: Optional[float] = Field(default=1200.0)
    kkanohm: Optional[float] = Field(default=0.4)
    ohm_threshsw: Optional[float] = Field(default=0.0)
    ohm_threshwd: Optional[float] = Field(default=0.0)
    ohm_coef: Optional[OHMCoefficients] = Field(default_factory=OHMCoefficients)
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

        # Get surface name
        surf_name = self.get_surface_name()

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = f"({surf_idx},)"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = np.nan
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # Get all properties of this class using introspection
        properties = [
            "sfr",
            "emis",
            "chanohm",
            "cpanohm",
            "kkanohm",
            "ohm_coef",
            "ohm_threshsw",
            "ohm_threshwd",
            "soildepth",
            "soilstorecap",
            "statelimit",
            "wetthresh",
            "sathydraulicconduct",
            "waterdist",
            "storedrainprm",
            "snowpacklimit",
            "thermal_layers",
            "irrfrac",
        ]
        # drop 'surface_type' and model-specific properties (e.g. model_xx)
        properties = [
            p for p in properties if p != "surface_type" and not p.startswith("model_")
        ]

        # Process each property
        dfs = [df_state]  # List to collect all DataFrames

        for property in properties:
            # Handle nested properties with their own to_df_state methods
            if property in [
                "waterdist",
                "storedrainprm",
                "ohm_coef",
                "lai",
            ]:
                nested_obj = getattr(self, property)
                if nested_obj is not None and hasattr(nested_obj, "to_df_state"):
                    nested_df = nested_obj.to_df_state(grid_id, surf_idx)
                    dfs.append(nested_df)
            elif property == "thermal_layers":
                nested_df = self.thermal_layers.to_df_state(
                    grid_id, surf_idx, surf_name
                )
                dfs.append(nested_df)
            elif property == "irrfrac":
                value = getattr(self, property)
                df_state.loc[grid_id, (f"{property}{surf_name}", "0")] = value
            elif property in ["sfr", "soilstorecap", "statelimit", "wetthresh"]:
                value = getattr(self, property)
                set_df_value(f"{property}_surf", value)
            else:
                value = getattr(self, property)
                set_df_value(property, value)
            # except Exception as e:
            #     print(f"Warning: Could not set property {property}: {str(e)}")
            #     continue

        # add dummy columns to conform to SUEWS convention
        list_cols = [
            "ohm_threshsw",
            "ohm_threshwd",
        ]
        for col in list_cols:
            df_state[(col, "(7,)")] = 0

        # Merge all DataFrames
        df_final = pd.concat(dfs, axis=1).sort_index(axis=1)
        return df_final

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "SurfaceProperties":
        """Reconstruct surface properties from DataFrame state format."""

        # Get surface index
        surf_idx = cls.get_surface_index()

        # Get surface name
        surf_name = cls.get_surface_name()

        # Get all properties of this class using introspection
        properties = [
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
            "waterdist",
            "storedrainprm",
            "snowpacklimit",
            "thermal_layers",
            "irrfrac",
        ]
        # drop 'surface_type' and model-specific properties (e.g. model_xx)
        properties = [
            p for p in properties if p != "surface_type" and not p.startswith("model_")
        ]

        # Process each property
        for property in properties:
            # Handle nested properties with their own from_df_state methods
            if property in [
                "waterdist",
                "storedrainprm",
                "ohm_coef",
                "lai",
            ]:
                nested_obj = getattr(cls, property)
                if nested_obj is not None and hasattr(nested_obj, "from_df_state"):
                    setattr(
                        cls, property, nested_obj.from_df_state(df, grid_id, surf_idx)
                    )
                continue
            elif property == "thermal_layers":
                setattr(
                    cls,
                    property,
                    cls.thermal_layers.from_df_state(df, grid_id, surf_idx, surf_name),
                )
            elif property == "irrfrac":
                setattr(
                    cls, property, df.loc[grid_id, (f"{property}{surf_name}", "0")]
                )
            elif property in ["sfr", "soilstorecap", "statelimit", "wetthresh"]:
                setattr(
                    cls,
                    property,
                    df.loc[grid_id, (f"{property}_surf", f"({surf_idx},)")],
                )
            else:
                setattr(cls, property, df.loc[grid_id, (property, f"({surf_idx},)")])

        return cls


class NonVegetatedSurfaceProperties(SurfaceProperties):
    alb: float = Field(ge=0, le=1, description="Surface albedo", default=0.1)

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert non-vegetated surface properties to DataFrame state format."""

        # Get base properties from parent
        df_base = super().to_df_state(grid_id)

        surf_idx = self.get_surface_index()

        if self.waterdist is not None:
            df_waterdist = self.waterdist.to_df_state(grid_id, surf_idx)
            df_base = pd.concat([df_base, df_waterdist], axis=1).sort_index(axis=1)

        for attr in ["alb"]:
            df_base.loc[grid_id, (attr, f"({surf_idx},)")] = getattr(self, attr)
            df_base = df_base.sort_index(axis=1)

        return df_base


class PavedProperties(NonVegetatedSurfaceProperties):
    _surface_type: Literal[SurfaceType.PAVED] = SurfaceType.PAVED
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.PAVED),
        description="Water distribution fractions for paved surfaces",
    )

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
                    "_surface_type",
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
                if nested_prop == "thermal_layers":
                    surf_name = self.get_surface_name()
                    nested_df = nested_obj.to_df_state(grid_id, surf_idx, surf_name)
                else:
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
    _facet_type: Literal["roof", "wall"] = PrivateAttr(default="roof")

    def to_df_state(
        self,
        grid_id: int,
        layer_idx: int,
        facet_type: Literal["roof", "wall"],
    ) -> pd.DataFrame:
        """Convert building layer parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            layer_idx: Layer index (0 or 1 for two layers)

        Returns:
            pd.DataFrame: DataFrame containing building layer parameters
        """
        df_state = init_df_state(grid_id)

        # Add basic parameters
        df_state[(f"alb_{facet_type}", f"({layer_idx},)")] = self.alb
        df_state[(f"emis_{facet_type}", f"({layer_idx},)")] = self.emis
        df_state[(f"statelimit_{facet_type}", f"({layer_idx},)")] = self.statelimit
        df_state[(f"soilstorecap_{facet_type}", f"({layer_idx},)")] = self.soilstorecap
        df_state[(f"wetthresh_{facet_type}", f"({layer_idx},)")] = self.wetthresh

        # Determine prefix based on layer type
        prefix = facet_type

        # Add layer-specific parameters
        if facet_type == "roof" and self.roof_albedo_dir_mult_fact is not None:
            df_state[(f"{prefix}_albedo_dir_mult_fact", f"(0, {layer_idx})")] = (
                self.roof_albedo_dir_mult_fact
            )
        elif facet_type == "wall" and self.wall_specular_frac is not None:
            df_state[(f"{prefix}_specular_frac", f"(0, {layer_idx})")] = (
                self.wall_specular_frac
            )

        # Add thermal layers
        df_thermal = self.thermal_layers.to_df_state(grid_id, layer_idx, facet_type)
        df_state = pd.concat([df_state, df_thermal], axis=1)

        return df_state

    @classmethod
    def from_df_state(
        cls,
        df: pd.DataFrame,
        grid_id: int,
        layer_idx: int,
        facet_type: Literal["roof", "wall"],
    ) -> "BuildingLayer":
        """Reconstruct BuildingLayer instance from DataFrame.

        Args:
            df: DataFrame containing building layer parameters.
            grid_id: Grid ID for the DataFrame index.
            layer_idx: Layer index (0 or 1 for two layers).
            facet_type: Facet type ("roof" or "wall").

        Returns:
            BuildingLayer: Reconstructed BuildingLayer instance.
        """
        # Prefix for the specific layer type
        prefix = facet_type

        # Extract scalar parameters
        params = {
            "alb": df.loc[grid_id, (f"alb_{prefix}", f"({layer_idx},)")],
            "emis": df.loc[grid_id, (f"emis_{prefix}", f"({layer_idx},)")],
            "statelimit": df.loc[grid_id, (f"statelimit_{prefix}", f"({layer_idx},)")],
            "soilstorecap": df.loc[
                grid_id, (f"soilstorecap_{prefix}", f"({layer_idx},)")
            ],
            "wetthresh": df.loc[grid_id, (f"wetthresh_{prefix}", f"({layer_idx},)")],
        }

        # Extract optional parameters
        if facet_type == "roof":
            params["roof_albedo_dir_mult_fact"] = df.loc[
                grid_id, (f"roof_albedo_dir_mult_fact", f"(0, {layer_idx})")]
            params["wall_specular_frac"] = None  # Explicitly set to None for clarity
        elif facet_type == "wall":
            params["wall_specular_frac"] = df.loc[
                grid_id, (f"wall_specular_frac", f"(0, {layer_idx})")]
            params["roof_albedo_dir_mult_fact"] = None  # Explicitly set to None for clarity

        # Extract ThermalLayers
        thermal_layers = ThermalLayers.from_df_state(df, grid_id, layer_idx, facet_type)

        # Add thermal_layers to params
        params["thermal_layers"] = thermal_layers

        # Return the reconstructed instance
        return cls(**params)


class RoofLayer(BuildingLayer):
    _facet_type: Literal["roof"] = "roof"


class WallLayer(BuildingLayer):
    _facet_type: Literal["wall"] = "wall"


class VerticalLayers(BaseModel):
    nlayer: int = Field(
        default=3, description="Number of vertical layers in the urban canopy"
    )
    height: List[float] = Field(
        default=[0.0, 10.0, 20.0, 30.0],
        description="Heights of layer boundaries in metres, length must be nlayer+1",
    )
    veg_frac: List[float] = Field(
        default=[0.0, 0.0, 0.0],
        description="Fraction of vegetation in each layer, length must be nlayer",
    )
    veg_scale: List[float] = Field(
        default=[1.0, 1.0, 1.0],
        description="Scaling factor for vegetation in each layer, length must be nlayer",
    )
    building_frac: List[float] = Field(
        default=[0.4, 0.3, 0.3],
        description="Fraction of buildings in each layer, must sum to 1.0, length must be nlayer",
    )
    building_scale: List[float] = Field(
        default=[1.0, 1.0, 1.0],
        description="Scaling factor for buildings in each layer, length must be nlayer",
    )
    roofs: List[RoofLayer] = Field(
        default_factory=lambda: [RoofLayer(), RoofLayer(), RoofLayer()],
        description="Properties for roof surfaces in each layer, length must be nlayer",
    )
    walls: List[WallLayer] = Field(
        default_factory=lambda: [WallLayer(), WallLayer(), WallLayer()],
        description="Properties for wall surfaces in each layer, length must be nlayer",
    )

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
        # This rule is not correct, we just need building_frac to be in range [0,1]
        # if not math.isclose(sum(self.building_frac), 1.0, rel_tol=1e-9):
        #    raise ValueError(
        #        f"Building fractions must sum to 1.0, got {sum(self.building_frac)}"
        #    )

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
            [self.roofs[i].to_df_state(grid_id, i, "roof") for i in range(self.nlayer)],
            axis=1,
        )
        df_walls = pd.concat(
            [self.walls[i].to_df_state(grid_id, i, "wall") for i in range(self.nlayer)],
            axis=1,
        )

        # Combine all DataFrames
        df_state = pd.concat([df_state, df_roofs, df_walls], axis=1)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "VerticalLayers":
        """Reconstruct VerticalLayers instance from DataFrame."""
        # Extract the number of layers
        nlayer = int(df.loc[grid_id, ("nlayer", "0")])

        # Extract heights for each layer boundary
        height = [df.loc[grid_id, ("height", f"({i},)")] for i in range(nlayer + 1)]

        # Extract vegetation and building parameters for each layer
        veg_frac = [df.loc[grid_id, ("veg_frac", f"({i},)")] for i in range(nlayer)]
        veg_scale = [df.loc[grid_id, ("veg_scale", f"({i},)")] for i in range(nlayer)]
        building_frac = [
            df.loc[grid_id, ("building_frac", f"({i},)")] for i in range(nlayer)]
        building_scale = [
            df.loc[grid_id, ("building_scale", f"({i},)")] for i in range(nlayer)]

        # Reconstruct roof and wall properties for each layer
        roofs = [RoofLayer.from_df_state(df, grid_id, i, "roof") for i in range(nlayer)]
        walls = [WallLayer.from_df_state(df, grid_id, i, "wall") for i in range(nlayer)]

        # Construct and return VerticalLayers instance
        return cls(
            nlayer=nlayer,
            height=height,
            veg_frac=veg_frac,
            veg_scale=veg_scale,
            building_frac=building_frac,
            building_scale=building_scale,
            roofs=roofs,
            walls=walls,
        )


class BldgsProperties(NonVegetatedSurfaceProperties):
    _surface_type: Literal[SurfaceType.BLDGS] = SurfaceType.BLDGS
    faibldg: float = Field(
        ge=0, default=0.3, description="Frontal area index of buildings"
    )
    bldgh: float = Field(ge=0, default=10.0, description="Building height")
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.BLDGS)
    )

    @model_validator(mode="after")
    def validate_rsl_zd_range(self) -> "BldgsProperties":
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
        df_state = super().to_df_state(grid_id).sort_index(axis=1)

        df_state.loc[grid_id, ("faibldg", "0")] = self.faibldg
        df_state = df_state.sort_index(axis=1)
        df_state.loc[grid_id, ("bldgh", "0")] = self.bldgh

        return df_state


class BsoilProperties(NonVegetatedSurfaceProperties):
    _surface_type: Literal[SurfaceType.BSOIL] = SurfaceType.BSOIL
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.BSOIL),
        description="Water distribution for bare soil",
    )

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert bare soil properties to DataFrame state format."""
        df_state = super().to_df_state(grid_id)
        # df_state.loc[grid_id, ("waterdist", "0")] = self.waterdist
        return df_state


class WaterProperties(NonVegetatedSurfaceProperties):
    _surface_type: Literal[SurfaceType.WATER] = SurfaceType.WATER
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

        list_attr = ["flowchange"]

        # Add all non-inherited properties
        df_state.loc[grid_id, ("flowchange", "0")] = self.flowchange

        return df_state


class ModelControl(BaseModel):
    tstep: int = Field(
        default=300, description="Time step in seconds for model calculations"
    )
    forcing_file: str = Field(
        default="forcing.txt", description="Path to meteorological forcing data file"
    )
    output_file: str = Field(
        default="output.txt", description="Path to model output file"
    )
    # daylightsaving_method: int
    diagnose: int = Field(
        default=0,
        description="Level of diagnostic output (0=none, 1=basic, 2=detailed)",
    )

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
    netradiationmethod: int = Field(
        default=3, description="Method used to calculate net radiation"
    )
    emissionsmethod: int = Field(
        default=2, description="Method used to calculate anthropogenic emissions"
    )
    storageheatmethod: int = Field(
        default=1, description="Method used to calculate storage heat flux"
    )
    ohmincqf: int = Field(
        default=0,
        description="Include anthropogenic heat in OHM calculations (1) or not (0)",
    )
    roughlenmommethod: int = Field(
        default=2, description="Method used to calculate momentum roughness length"
    )
    roughlenheatmethod: int = Field(
        default=2, description="Method used to calculate heat roughness length"
    )
    stabilitymethod: int = Field(
        default=2, description="Method used for atmospheric stability calculation"
    )
    smdmethod: int = Field(
        default=1, description="Method used to calculate soil moisture deficit"
    )
    waterusemethod: int = Field(
        default=1, description="Method used to calculate water use"
    )
    diagmethod: int = Field(default=1, description="Method used for model diagnostics")
    faimethod: int = Field(
        default=1, description="Method used to calculate frontal area index"
    )
    localclimatemethod: int = Field(
        default=0, description="Method used for local climate zone calculations"
    )
    snowuse: int = Field(
        default=0, description="Include snow calculations (1) or not (0)"
    )

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
    rainmaxres: float = Field(ge=0, le=20, default=0.25)
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
    air_ext_lw: float = Field(
        default=0.0, description="Air extinction coefficient for longwave radiation"
    )
    air_ext_sw: float = Field(
        default=0.0, description="Air extinction coefficient for shortwave radiation"
    )
    air_ssa_lw: float = Field(
        default=0.5, description="Air single scattering albedo for longwave radiation"
    )
    air_ssa_sw: float = Field(
        default=0.5, description="Air single scattering albedo for shortwave radiation"
    )
    ground_albedo_dir_mult_fact: float = Field(
        default=1.0, description="Multiplication factor for direct ground albedo"
    )
    n_stream_lw_urban: int = Field(
        default=2, description="Number of streams for longwave radiation in urban areas"
    )
    n_stream_sw_urban: int = Field(
        default=2,
        description="Number of streams for shortwave radiation in urban areas",
    )
    n_vegetation_region_urban: int = Field(
        default=1, description="Number of vegetation regions in urban areas"
    )
    sw_dn_direct_frac: float = Field(
        default=0.5,
        description="Fraction of downward shortwave radiation that is direct",
    )
    use_sw_direct_albedo: float = Field(
        default=1.0, description="Flag to use direct albedo for shortwave radiation"
    )
    veg_contact_fraction_const: float = Field(
        default=0.5, description="Constant vegetation contact fraction"
    )
    veg_fsd_const: float = Field(
        default=0.5, description="Constant vegetation fractional standard deviation"
    )
    veg_ssa_lw: float = Field(
        default=0.5,
        description="Vegetation single scattering albedo for longwave radiation",
    )
    veg_ssa_sw: float = Field(
        default=0.5,
        description="Vegetation single scattering albedo for shortwave radiation",
    )

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
    working_day: float = Field(default=1.0)
    holiday: float = Field(default=0.0)

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
    monday: float = 0.0
    tuesday: float = 0.0
    wednesday: float = 0.0
    thursday: float = 0.0
    friday: float = 0.0
    saturday: float = 0.0
    sunday: float = 0.0

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

    @classmethod
    def __init_default_values__(cls) -> Dict[str, Dict[str, float]]:
        """Generate default values for hourly profiles.

        Returns:
            Dict containing default working_day and holiday profiles with uniform distribution
        """
        # Create uniform distribution (1/24) for each hour
        uniform_value = 1.0 / 24.0

        # Generate hour keys 1-24 with uniform values
        hourly_values = {str(hour): uniform_value for hour in range(1, 25)}

        return {"working_day": hourly_values.copy(), "holiday": hourly_values.copy()}

    def __init__(self, **data):
        # If no values provided, use defaults
        if not data:
            defaults = self.__init_default_values__()
            data = defaults
        super().__init__(**data)

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
            df_state[(param_name, f"({int(hour)-1}, 0)")] = value

        # Set holiday/weekend values (index 1)
        for hour, value in self.holiday.items():
            df_state[(param_name, f"({int(hour)-1}, 1)")] = value

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
            str(hour + 1): df.loc[grid_id, (param_name, f"({hour}, 0)")]
            for hour in range(24)
        }

        # Extract holiday values (index 1)
        holiday = {
            str(hour + 1): df.loc[grid_id, (param_name, f"({hour}, 1)")]
            for hour in range(24)
        }

        # Create an instance of HourlyProfile
        return cls(working_day=working_day, holiday=holiday)


class IrrigationParams(BaseModel):
    h_maintain: float = Field(
        default=0.5, description="Soil moisture threshold for irrigation"
    )
    faut: float = Field(default=0.0, description="Fraction of automatic irrigation")
    ie_start: float = Field(default=0.0, description="Start time of irrigation (hour)")
    ie_end: float = Field(default=0.0, description="End time of irrigation (hour)")
    internalwateruse_h: float = Field(
        default=0.0, description="Internal water use per hour"
    )
    daywatper: WeeklyProfile = Field(default_factory=WeeklyProfile)
    daywat: WeeklyProfile = Field(default_factory=WeeklyProfile)
    wuprofa_24hr: HourlyProfile = Field(default_factory=HourlyProfile)
    wuprofm_24hr: HourlyProfile = Field(default_factory=HourlyProfile)

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """
        Convert irrigation parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing irrigation parameters
        """

        df_state = init_df_state(grid_id)

        df_state.loc[grid_id, ("h_maintain", "0")] = self.h_maintain
        df_state.loc[grid_id, ("faut", "0")] = self.faut
        df_state.loc[grid_id, ("ie_start", "0")] = self.ie_start
        df_state.loc[grid_id, ("ie_end", "0")] = self.ie_end
        df_state.loc[grid_id, ("internalwateruse_h", "0")] = self.internalwateruse_h

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
        h_maintain = df.loc[grid_id, ("h_maintain", "0")]
        faut = df.loc[grid_id, ("faut", "0")]
        ie_start = df.loc[grid_id, ("ie_start", "0")]
        ie_end = df.loc[grid_id, ("ie_end", "0")]
        internalwateruse_h = df.loc[grid_id, ("internalwateruse_h", "0")]

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
    qf0_beu: DayProfile = Field(
        description="Base anthropogenic heat flux for buildings, equipment and urban metabolism",
        default_factory=DayProfile,
    )
    qf_a: DayProfile = Field(
        description="Coefficient a for anthropogenic heat flux calculation",
        default_factory=DayProfile,
    )
    qf_b: DayProfile = Field(
        description="Coefficient b for anthropogenic heat flux calculation",
        default_factory=DayProfile,
    )
    qf_c: DayProfile = Field(
        description="Coefficient c for anthropogenic heat flux calculation",
        default_factory=DayProfile,
    )
    baset_cooling: DayProfile = Field(
        description="Base temperature for cooling degree days",
        default_factory=DayProfile,
    )
    baset_heating: DayProfile = Field(
        description="Base temperature for heating degree days",
        default_factory=DayProfile,
    )
    ah_min: DayProfile = Field(
        description="Minimum anthropogenic heat flux", default_factory=DayProfile
    )
    ah_slope_cooling: DayProfile = Field(
        description="Slope of anthropogenic heat vs cooling degree days",
        default_factory=DayProfile,
    )
    ah_slope_heating: DayProfile = Field(
        description="Slope of anthropogenic heat vs heating degree days",
        default_factory=DayProfile,
    )
    ahprof_24hr: HourlyProfile = Field(
        description="24-hour profile of anthropogenic heat flux",
        default_factory=HourlyProfile,
    )
    popdensdaytime: DayProfile = Field(
        description="Daytime population density", default_factory=DayProfile
    )
    popdensnighttime: float = Field(
        default=10.0, description="Nighttime population density"
    )
    popprof_24hr: HourlyProfile = Field(
        description="24-hour profile of population density",
        default_factory=HourlyProfile,
    )

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

        df_state.loc[grid_id, ("popdensnighttime", "0")] = self.popdensnighttime

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
        popdensnighttime = df.loc[grid_id, ("popdensnighttime", "0")]

        # Construct and return AnthropogenicHeat instance
        return cls(
            **day_profiles,
            **hourly_profiles,
            popdensnighttime=popdensnighttime,
        )


class CO2Params(BaseModel):
    co2pointsource: float = Field(
        default=0.0, description="CO2 point source emission factor"
    )
    ef_umolco2perj: float = Field(
        default=0.0, description="CO2 emission factor per unit of fuel"
    )
    enef_v_jkm: float = Field(
        default=0.0, description="CO2 emission factor per unit of vehicle distance"
    )
    fcef_v_kgkm: DayProfile = Field(
        description="Fuel consumption efficiency for vehicles",
        default_factory=DayProfile,
    )
    frfossilfuel_heat: float = Field(
        default=0.0, description="Fraction of fossil fuel heat"
    )
    frfossilfuel_nonheat: float = Field(
        default=0.0, description="Fraction of fossil fuel non-heat"
    )
    maxfcmetab: float = Field(
        default=0.0, description="Maximum fuel consumption metabolic rate"
    )
    maxqfmetab: float = Field(
        default=0.0, description="Maximum heat production metabolic rate"
    )
    minfcmetab: float = Field(
        default=0.0, description="Minimum fuel consumption metabolic rate"
    )
    minqfmetab: float = Field(
        default=0.0, description="Minimum heat production metabolic rate"
    )
    trafficrate: DayProfile = Field(
        description="Traffic rate", default_factory=DayProfile
    )
    trafficunits: float = Field(default=0.0, description="Traffic units")
    traffprof_24hr: HourlyProfile = Field(
        description="24-hour profile of traffic rate", default_factory=HourlyProfile
    )
    humactivity_24hr: HourlyProfile = Field(
        description="24-hour profile of human activity", default_factory=HourlyProfile
    )

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
            df_state.loc[grid_id, (param_name, "0")] = value

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
            "co2pointsource": df.loc[grid_id, ("co2pointsource", "0")],
            "ef_umolco2perj": df.loc[grid_id, ("ef_umolco2perj", "0")],
            "enef_v_jkm": df.loc[grid_id, ("enef_v_jkm", "0")],
            "frfossilfuel_heat": df.loc[grid_id, ("frfossilfuel_heat", "0")],
            "frfossilfuel_nonheat": df.loc[grid_id, ("frfossilfuel_nonheat", "0")],
            "maxfcmetab": df.loc[grid_id, ("maxfcmetab", "0")],
            "maxqfmetab": df.loc[grid_id, ("maxqfmetab", "0")],
            "minfcmetab": df.loc[grid_id, ("minfcmetab", "0")],
            "minqfmetab": df.loc[grid_id, ("minqfmetab", "0")],
            "trafficunits": df.loc[grid_id, ("trafficunits", "0")],
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
    startdls: float = Field(
        default=0.0, description="Start of daylight savings time in decimal day of year"
    )
    enddls: float = Field(
        default=0.0, description="End of daylight savings time in decimal day of year"
    )
    heat: AnthropogenicHeat = Field(
        description="Anthropogenic heat emission parameters",
        default_factory=AnthropogenicHeat,
    )
    co2: CO2Params = Field(
        description="CO2 emission parameters", default_factory=CO2Params
    )

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
        df_state.loc[grid_id, ("startdls", "0")] = self.startdls
        df_state.loc[grid_id, ("enddls", "0")] = self.enddls

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
        startdls = df.loc[grid_id, ("startdls", "0")]
        enddls = df.loc[grid_id, ("enddls", "0")]

        # Reconstruct heat parameters
        heat = AnthropogenicHeat.from_df_state(df, grid_id)

        # Reconstruct CO2 parameters
        co2 = CO2Params.from_df_state(df, grid_id)

        return cls(startdls=startdls, enddls=enddls, heat=heat, co2=co2)


class Conductance(BaseModel):
    g_max: float = Field(default=40.0, description="Maximum conductance")
    g_k: float = Field(
        default=0.6,
        description="Conductance parameter related to incoming solar radiation",
    )
    g_q_base: float = Field(
        default=0.03,
        description="Base value for conductance parameter related to vapor pressure deficit",
    )
    g_q_shape: float = Field(
        default=0.9,
        description="Shape parameter for conductance related to vapor pressure deficit",
    )
    g_t: float = Field(
        default=30.0, description="Conductance parameter related to air temperature"
    )
    g_sm: float = Field(
        default=0.5, description="Conductance parameter related to soil moisture"
    )
    kmax: float = Field(
        default=1200.0, description="Maximum incoming shortwave radiation"
    )
    gsmodel: int = Field(default=1, description="Stomatal conductance model selection")
    s1: float = Field(default=0.2, description="Soil moisture threshold parameter")
    s2: float = Field(default=0.5, description="Soil moisture threshold parameter")
    tl: float = Field(default=0.0, description="Air temperature threshold parameter")
    th: float = Field(default=50.0, description="Air temperature threshold parameter")

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
            df_state.loc[grid_id, (param_name, "0")] = value

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
            "gsmodel": int(df.loc[grid_id, ("gsmodel", "0")]),
            "s1": df.loc[grid_id, ("s1", "0")],
            "s2": df.loc[grid_id, ("s2", "0")],
            "tl": df.loc[grid_id, ("tl", "0")],
            "th": df.loc[grid_id, ("th", "0")],
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
            df.loc[grid_id, ("laipower", f"(0, {veg_idx})")],
            df.loc[grid_id, ("laipower", f"(1, {veg_idx})")],
            df.loc[grid_id, ("laipower", f"(2, {veg_idx})")],
            df.loc[grid_id, ("laipower", f"(3, {veg_idx})")],
        ]

        # Return the instance with coefficients
        return cls(
            growth_lai=coefficients[0],
            growth_gdd=coefficients[1],
            senescence_lai=coefficients[2],
            senescence_sdd=coefficients[3],
        )


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

        # Extract LAI power coefficients
        laipower = LAIPowerCoefficients.from_df_state(df, grid_id, veg_idx)

        return cls(**lai_params, laipower=laipower)


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
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = np.nan
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # add ordinary float properties
        for attr in [
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
            set_df_value(attr, f"({surf_idx-2},)", getattr(self, attr))

        df_lai = self.lai.to_df_state(grid_id, surf_idx)
        df_state = pd.concat([df_state, df_lai], axis=1).sort_index(axis=1)

        return df_state


class DectrProperties(VegetatedSurfaceProperties):
    faidectree: float = Field(
        default=0.1, description="Frontal area index of deciduous trees"
    )
    dectreeh: float = Field(default=15.0, description="Deciduous tree height")
    pormin_dec: float = Field(default=0.2, description="Minimum porosity")
    pormax_dec: float = Field(default=0.6, description="Maximum porosity")
    capmax_dec: float = Field(default=100.0, description="Maximum capacity")
    capmin_dec: float = Field(default=10.0, description="Minimum capacity")
    _surface_type: Literal[SurfaceType.DECTR] = SurfaceType.DECTR
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.DECTR),
        description="Water distribution for deciduous trees",
    )

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
            df_state.loc[grid_id, (attr, "0")] = getattr(self, attr)

        # specific properties
        df_state.loc[grid_id, ("albmin_dectr", "0")] = self.alb_min
        df_state.loc[grid_id, ("albmax_dectr", "0")] = self.alb_max

        return df_state


class EvetrProperties(VegetatedSurfaceProperties):
    faievetree: float = Field(
        default=0.1, description="Frontal area index of evergreen trees"
    )
    evetreeh: float = Field(default=15.0, description="Evergreen tree height")
    _surface_type: Literal[SurfaceType.EVETR] = SurfaceType.EVETR
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.EVETR),
        description="Water distribution for evergreen trees",
    )

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
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # Add all non-inherited properties
        list_properties = ["faievetree", "evetreeh"]
        for attr in list_properties:
            df_state.loc[grid_id, (attr, "0")] = getattr(self, attr)

        # specific properties
        df_state.loc[grid_id, ("albmin_evetr", "0")] = self.alb_min
        df_state.loc[grid_id, ("albmax_evetr", "0")] = self.alb_max

        return df_state


class GrassProperties(VegetatedSurfaceProperties):
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
        df_state[("albmin_grass", "0")] = self.alb_min
        df_state[("albmax_grass", "0")] = self.alb_max

        return df_state


class SnowParams(BaseModel):
    crwmax: float = Field(default=0.1, description="Maximum water capacity of snow")
    crwmin: float = Field(default=0.05, description="Minimum water capacity of snow")
    narp_emis_snow: float = Field(default=0.99, description="Snow surface emissivity")
    preciplimit: float = Field(
        default=2.2, description="Limit for snow vs rain precipitation"
    )
    preciplimitalb: float = Field(
        default=0.1, description="Precipitation limit for albedo aging"
    )
    snowalbmax: float = Field(default=0.85, description="Maximum snow albedo")
    snowalbmin: float = Field(default=0.4, description="Minimum snow albedo")
    snowdensmin: float = Field(
        default=100.0, description="Minimum snow density (kg m-3)"
    )
    snowdensmax: float = Field(
        default=400.0, description="Maximum snow density (kg m-3)"
    )
    snowlimbldg: float = Field(default=0.1, description="Snow limit on buildings")
    snowlimpaved: float = Field(default=0.1, description="Snow limit on paved surfaces")
    snowprof_24hr: HourlyProfile = Field(
        default_factory=HourlyProfile, description="24-hour snow profile"
    )
    tau_a: float = Field(default=0.018, description="Aging constant for cold snow")
    tau_f: float = Field(default=0.11, description="Aging constant for melting snow")
    tau_r: float = Field(default=0.05, description="Aging constant for refreezing snow")
    tempmeltfact: float = Field(default=0.12, description="Temperature melt factor")
    radmeltfact: float = Field(default=0.0016, description="Radiation melt factor")

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
            df_state.loc[grid_id, (param_name, "0")] = value

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
    dectr: DectrProperties = Field(
        default_factory=DectrProperties,
        description="Properties for deciduous trees and vegetation",
    )
    evetr: EvetrProperties = Field(
        default_factory=EvetrProperties,
        description="Properties for evergreen trees and vegetation",
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
            "dectr": DectrProperties.from_df_state(df, grid_id),
            "evetr": EvetrProperties.from_df_state(df, grid_id),
            "grass": GrassProperties.from_df_state(df, grid_id),
            "bsoil": BsoilProperties.from_df_state(df, grid_id),
            "water": WaterProperties.from_df_state(df, grid_id)
        }

        # Return reconstructed instance
        return cls(**params)


class SiteProperties(BaseModel):
    lat: float = Field(
        ge=-90, le=90, description="Latitude of the site in degrees", default=51.5
    )
    lng: float = Field(
        ge=-180, le=180, description="Longitude of the site in degrees", default=-0.13
    )
    alt: float = Field(
        gt=0, description="Altitude of the site in metres above sea level", default=40.0
    )
    timezone: int = Field(
        ge=-12, le=12, description="Time zone offset from UTC in hours", default=0
    )
    surfacearea: float = Field(
        gt=0,
        description="Total surface area of the site in square metres",
        default=10000.0,
    )
    z: float = Field(gt=0, description="Measurement height in metres", default=10.0)
    z0m_in: float = Field(
        gt=0, description="Momentum roughness length in metres", default=1.0
    )
    zdm_in: float = Field(
        gt=0, description="Zero-plane displacement height in metres", default=5.0
    )
    pipecapacity: float = Field(
        gt=0, description="Maximum capacity of drainage pipes in mm/hr", default=100.0
    )
    runofftowater: float = Field(
        ge=0,
        le=1,
        description="Fraction of excess water going to water bodies",
        default=0.0,
    )
    narp_trans_site: float = Field(
        description="Site-specific NARP transmission coefficient", default=0.2
    )
    lumps: LUMPSParams = Field(
        default_factory=LUMPSParams,
        description="Parameters for Local-scale Urban Meteorological Parameterization Scheme",
    )
    spartacus: SPARTACUSParams = Field(
        default_factory=SPARTACUSParams,
        description="Parameters for Solar Parametrizations for Radiative Transfer through Urban Canopy Scheme",
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
        ]:
            df_state.loc[grid_id, (f"{var}", "0")] = getattr(self, var)

        # complex attributes
        df_lumps = self.lumps.to_df_state(grid_id)
        df_spartacus = self.spartacus.to_df_state(grid_id)
        df_conductance = self.conductance.to_df_state(grid_id)
        df_irrigation = self.irrigation.to_df_state(grid_id)
        df_anthropogenic_emissions = self.anthropogenic_emissions.to_df_state(grid_id)
        df_snow = self.snow.to_df_state(grid_id)
        df_land_cover = self.land_cover.to_df_state(grid_id)
        df_vertical_layers = self.vertical_layers.to_df_state(grid_id)

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
        ]:
            params[var] = df.loc[grid_id, (var, "0")]

        # Extract complex attributes
        params["lumps"] = LUMPSParams.from_df_state(df, grid_id)
        params["spartacus"] = SPARTACUSParams.from_df_state(df, grid_id)
        params["conductance"] = Conductance.from_df_state(df, grid_id)
        params["irrigation"] = IrrigationParams.from_df_state(df, grid_id)
        params["anthropogenic_emissions"] = AnthropogenicEmissions.from_df_state(df, grid_id)
        params["snow"] = SnowParams.from_df_state(df, grid_id)
        params["land_cover"] = LandCover.from_df_state(df, grid_id)
        params["vertical_layers"] = VerticalLayers.from_df_state(df, grid_id)

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


class Model(BaseModel):
    control: ModelControl = Field(
        default_factory=ModelControl,
        description="Model control parameters including timestep, output options, etc.",
    )
    physics: ModelPhysics = Field(
        default_factory=ModelPhysics,
        description="Model physics parameters including surface properties, coefficients, etc.",
    )

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert model to DataFrame state format"""
        df_state = init_df_state(grid_id)
        df_control = self.control.to_df_state(grid_id)
        df_physics = self.physics.to_df_state(grid_id)
        df_state = pd.concat([df_state, df_control, df_physics], axis=1)
        return df_state


class SUEWSConfig(BaseModel):
    name: str = Field(
        default="sample config", description="Name of the SUEWS configuration"
    )
    description: str = Field(
        default="this is a sample config for testing purposes ONLY - values are not realistic",
        description="Description of this SUEWS configuration",
    )
    model: Model = Field(
        default_factory=Model,
        description="Model control and physics parameters",
    )
    site: List[Site] = Field(
        default=[Site()],
        description="List of sites to simulate",
        min_items=1,
    )

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
        list_df_site = []
        for grid_id in range(len(self.site)):
            df_site = self.site[grid_id].to_df_state(grid_id)
            df_model = self.model.to_df_state(grid_id)
            df_site = pd.concat([df_site, df_model], axis=1)
            list_df_site.append(df_site)

        df = pd.concat(list_df_site, axis=1)
        # remove duplicate columns
        df = df.loc[:, ~df.columns.duplicated()]
        return df

    @classmethod
    def from_df_state(cls, df: pd.DataFrame) -> "SUEWSConfig":
        """Create config from DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing SUEWS configuration state.

        Returns:
            SUEWSConfig: Instance of SUEWSConfig reconstructed from DataFrame.
        """
        # Initialize with default values
        config = cls()

        # Get grid IDs from DataFrame index
        grid_ids = df.index.tolist()

        # Create list of sites
        sites = []
        for grid_id in grid_ids:
            # Create site instance
            site = Site(gridiv=grid_id)

            # Set site properties
            site_properties = SiteProperties.from_df_state(df, grid_id)
            site.properties = site_properties

            # Set initial states
            initial_states = InitialStates.from_df_state(df, grid_id)
            site.initial_states = initial_states

            sites.append(site)

        # Update config with reconstructed data
        config.site = sites

        # Reconstruct model
        model = Model()
        for grid_id in grid_ids:
            # Set model control
            model_control = ModelControl.from_df_state(df, grid_id)
            model.control = model_control

            # Set model physics
            model_physics = ModelPhysics.from_df_state(df, grid_id)
            model.physics = model_physics
            break  # Only need one as model is shared across sites

        config.model = model

        return config
