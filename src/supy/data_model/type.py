from typing import TypeVar, Optional, Generic
from pydantic import BaseModel
import numpy as np
import pandas as pd
from enum import Enum

class SurfaceType(str, Enum):
    PAVED = "paved"
    BLDGS = "bldgs"
    EVETR = "evetr"
    DECTR = "dectr"
    GRASS = "grass"
    BSOIL = "bsoil"
    WATER = "water"


T = TypeVar("T")


class Reference(BaseModel):
    desc: Optional[str] = None
    ID: Optional[str] = None
    DOI: Optional[str] = None


class ValueWithDOI(BaseModel, Generic[T]):
    """A class that wraps a value with an optional DOI reference.

    This class allows storing a value along with its reference information (e.g. DOI).
    It handles numeric type conversion and implements comparison operators.

    Attributes:
        value (T): The wrapped value of generic type T
        ref (Optional[Reference]): Optional reference information for the value
    """

    value: T
    ref: Optional[Reference] = None

    def __init__(self, value: T, ref: Optional[Reference] = None):
        # Convert numpy numeric types to Python native types
        if isinstance(value, (np.float64, np.float32)):
            value = float(value)
        elif isinstance(value, (np.int64, np.int32)):
            value = int(value)
        super().__init__(value=value, ref=ref)

    def __str__(self):
        """String representation showing just the value"""
        return f"{self.value}"

    def __repr__(self):
        """Representation showing just the value"""
        return f"{self.value}"

    def __eq__(self, other):
        """Equal comparison operator"""
        if isinstance(other, ValueWithDOI):
            return self.value == other.value
        return self.value == other

    def __lt__(self, other):
        """Less than comparison operator"""
        if isinstance(other, ValueWithDOI):
            return self.value < other.value
        return self.value < other

    def __le__(self, other):
        """Less than or equal comparison operator"""
        if isinstance(other, ValueWithDOI):
            return self.value <= other.value
        return self.value <= other

    def __gt__(self, other):
        """Greater than comparison operator"""
        if isinstance(other, ValueWithDOI):
            return self.value > other.value
        return self.value > other

    def __ge__(self, other):
        """Greater than or equal comparison operator"""
        if isinstance(other, ValueWithDOI):
            return self.value >= other.value
        return self.value >= other

    def __ne__(self, other):
        """Not equal comparison operator"""
        if isinstance(other, ValueWithDOI):
            return self.value != other.value
        return self.value != other


def init_df_state(grid_id: int) -> pd.DataFrame:
    idx = pd.Index([grid_id], name="grid")
    col = pd.MultiIndex.from_tuples([("gridiv", "0")], names=["var", "ind_dim"])
    df_state = pd.DataFrame(index=idx, columns=col)
    df_state.loc[grid_id, ("gridiv", "0")] = grid_id
    return df_state