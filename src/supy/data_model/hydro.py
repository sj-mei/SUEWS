from pydantic import BaseModel, Field, PrivateAttr
from typing import Optional
import pandas as pd
from .type import ValueWithDOI, Reference, SurfaceType
import math


class WaterDistribution(BaseModel):
    # Optional fields for all possible distributions
    to_paved: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_bldgs: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_dectr: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_evetr: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_grass: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_bsoil: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_water: Optional[ValueWithDOI[float]] = Field(None, ge=0, le=1)
    to_runoff: Optional[ValueWithDOI[float]] = Field(
        None, ge=0, le=1
    )  # For paved/bldgs
    to_soilstore: Optional[ValueWithDOI[float]] = Field(
        None, ge=0, le=1
    )  # For vegetated surfaces
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
            SurfaceType.EVETR: {
                "to_paved": ValueWithDOI(0.1),
                "to_bldgs": ValueWithDOI(0.1),
                "to_dectr": ValueWithDOI(0.1),
                "to_grass": ValueWithDOI(0.1),
                "to_bsoil": ValueWithDOI(0.1),
                "to_water": ValueWithDOI(0.1),
                "to_soilstore": ValueWithDOI(0.4),
            },
            SurfaceType.DECTR: {
                "to_paved": ValueWithDOI(0.1),
                "to_bldgs": ValueWithDOI(0.1),
                "to_evetr": ValueWithDOI(0.1),
                "to_grass": ValueWithDOI(0.1),
                "to_bsoil": ValueWithDOI(0.1),
                "to_water": ValueWithDOI(0.1),
                "to_soilstore": ValueWithDOI(0.4),
            },
            SurfaceType.GRASS: {
                "to_paved": ValueWithDOI(0.1),
                "to_bldgs": ValueWithDOI(0.1),
                "to_dectr": ValueWithDOI(0.1),
                "to_evetr": ValueWithDOI(0.1),
                "to_bsoil": ValueWithDOI(0.1),
                "to_water": ValueWithDOI(0.1),
                "to_soilstore": ValueWithDOI(0.4),
            },
            SurfaceType.BSOIL: {
                "to_paved": ValueWithDOI(0.1),
                "to_bldgs": ValueWithDOI(0.1),
                "to_dectr": ValueWithDOI(0.1),
                "to_evetr": ValueWithDOI(0.1),
                "to_grass": ValueWithDOI(0.1),
                "to_water": ValueWithDOI(0.1),
                "to_soilstore": ValueWithDOI(0.4),
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
        total = sum(
            value.value if isinstance(value, ValueWithDOI) else value
            for value in values
        )
        # if not np.isclose(total, 1.0, rtol=1e-5):
        if not math.isclose(total, 1.0, rel_tol=1e-5):
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

        # Convert ValueWithDOI to float
        values = [
            value.value if isinstance(value, ValueWithDOI) else value
            for value in values
        ]

        # Create DataFrame with single row
        df = pd.DataFrame(
            index=pd.Index([grid_id], name="grid"),
            columns=columns,
            data=[values],
            dtype=float,
        )

        return df

    @classmethod
    def from_df_state(
        cls, df: pd.DataFrame, grid_id: int, surf_idx: int
    ) -> "WaterDistribution":
        """
        Reconstruct WaterDistribution from a DataFrame state format.

        Args:
            df (pd.DataFrame): DataFrame containing water distribution parameters.
            grid_id (int): Grid ID for the DataFrame index.
            surf_idx (int): Surface index for identifying columns.

        Returns:
            WaterDistribution: Instance of WaterDistribution.
        """
        dict_surface_type = {
            0: SurfaceType.PAVED,
            1: SurfaceType.BLDGS,
            2: SurfaceType.EVETR,
            3: SurfaceType.DECTR,
            4: SurfaceType.GRASS,
            5: SurfaceType.BSOIL,
            6: SurfaceType.WATER,
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
            # "to_soilstore": 7,
            # "to_runoff": 8,
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
                    "drain_eq",
                    "drain_coef_1",
                    "drain_coef_2",
                    "store_max",
                    "store_cap",
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
                "drain_eq",
                "drain_coef_1",
                "drain_coef_2",
                "store_max",
                "store_cap",
            ]
        ):
            df.loc[grid_id, ("storedrainprm", f"({i}, {surf_idx})")] = getattr(
                self, var
            ).value

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
            "drain_eq": 1,
            "drain_coef_1": 2,
            "drain_coef_2": 3,
            "store_max": 4,
            "store_cap": 5,
        }

        # Extract the values from the DataFrame
        params = {
            param: df.loc[grid_id, ("storedrainprm", f"({idx}, {surf_idx})")]
            for param, idx in param_map.items()
        }

        # Conver params to ValueWithDOI
        params = {key: ValueWithDOI(value) for key, value in params.items()}

        # Create an instance using the extracted parameters
        return cls(**params)
