"""Base classes and types for surface models."""

from enum import Enum
from typing import Dict, List, Optional, Union, Literal
from pydantic import BaseModel, Field, model_validator, PrivateAttr
import pandas as pd

from ..base import ValueWithDOI, Reference


def init_df_state(grid_id: int) -> pd.DataFrame:
    """Initialize an empty DataFrame with proper index and column structure."""
    idx = pd.Index([grid_id], name="grid")
    col = pd.MultiIndex.from_tuples([("gridiv", "0")], names=["var", "ind_dim"])
    df_state = pd.DataFrame(index=idx, columns=col)
    df_state.loc[grid_id, ("gridiv", "0")] = grid_id
    return df_state


class SurfaceType(str, Enum):
    """Enumeration of surface types in SUEWS."""
    PAVED = "paved"
    BLDGS = "bldgs"
    EVETR = "evetr"
    DECTR = "dectr"
    GRASS = "grass"
    BSOIL = "bsoil"
    WATER = "water"
