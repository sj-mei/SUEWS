from pydantic import BaseModel, Field
from typing import Optional
import pandas as pd
from .type import ValueWithDOI, Reference
from .profile import HourlyProfile, WeeklyProfile, DayProfile
from .type import init_df_state


class IrrigationParams(
    BaseModel
):  # TODO: May need to add ValueWithDOI to the profiles here
    h_maintain: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.5), description="Soil moisture threshold for irrigation"
    )
    faut: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Fraction of automatic irrigation"
    )
    ie_start: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Start time of irrigation (hour)"
    )
    ie_end: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="End time of irrigation (hour)"
    )
    internalwateruse_h: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Internal water use per hour"
    )
    daywatper: WeeklyProfile = Field(default_factory=WeeklyProfile)
    daywat: WeeklyProfile = Field(default_factory=WeeklyProfile)
    wuprofa_24hr: HourlyProfile = Field(default_factory=HourlyProfile)
    wuprofm_24hr: HourlyProfile = Field(default_factory=HourlyProfile)

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """
        Convert irrigation parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index

        Returns:
            pd.DataFrame: DataFrame containing irrigation parameters
        """

        df_state = init_df_state(grid_id)

        df_state.loc[grid_id, ("h_maintain", "0")] = self.h_maintain.value
        df_state.loc[grid_id, ("faut", "0")] = self.faut.value
        df_state.loc[grid_id, ("ie_start", "0")] = self.ie_start.value
        df_state.loc[grid_id, ("ie_end", "0")] = self.ie_end.value
        df_state.loc[grid_id, ("internalwateruse_h", "0")] = (
            self.internalwateruse_h.value
        )

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

        # Conver to ValueWithDOI
        h_maintain = ValueWithDOI(h_maintain)
        faut = ValueWithDOI(faut)
        ie_start = ValueWithDOI(ie_start)
        ie_end = ValueWithDOI(ie_end)
        internalwateruse_h = ValueWithDOI(internalwateruse_h)

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


class AnthropogenicHeat(
    BaseModel
):  # TODO: May need to add the ValueWithDOI to the profiles here
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

    ref: Optional[Reference] = None

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


class CO2Params(
    BaseModel
):  # TODO: May need to add the ValueWithDOI to the profiles here
    co2pointsource: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="CO2 point source emission factor"
    )
    ef_umolco2perj: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="CO2 emission factor per unit of fuel"
    )
    enef_v_jkm: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0),
        description="CO2 emission factor per unit of vehicle distance",
    )
    fcef_v_kgkm: DayProfile = Field(
        description="Fuel consumption efficiency for vehicles",
        default_factory=DayProfile,
    )
    frfossilfuel_heat: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Fraction of fossil fuel heat"
    )
    frfossilfuel_nonheat: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Fraction of fossil fuel non-heat"
    )
    maxfcmetab: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Maximum fuel consumption metabolic rate"
    )
    maxqfmetab: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Maximum heat production metabolic rate"
    )
    minfcmetab: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Minimum fuel consumption metabolic rate"
    )
    minqfmetab: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Minimum heat production metabolic rate"
    )
    trafficrate: DayProfile = Field(
        description="Traffic rate", default_factory=DayProfile
    )
    trafficunits: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Traffic units"
    )
    traffprof_24hr: HourlyProfile = Field(
        description="24-hour profile of traffic rate", default_factory=HourlyProfile
    )
    humactivity_24hr: HourlyProfile = Field(
        description="24-hour profile of human activity", default_factory=HourlyProfile
    )

    ref: Optional[Reference] = None

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
            "co2pointsource": self.co2pointsource.value,
            "ef_umolco2perj": self.ef_umolco2perj.value,
            "enef_v_jkm": self.enef_v_jkm.value,
            "frfossilfuel_heat": self.frfossilfuel_heat.value,
            "frfossilfuel_nonheat": self.frfossilfuel_nonheat.value,
            "maxfcmetab": self.maxfcmetab.value,
            "maxqfmetab": self.maxqfmetab.value,
            "minfcmetab": self.minfcmetab.value,
            "minqfmetab": self.minqfmetab.value,
            "trafficunits": self.trafficunits.value,
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

        # Convert scalar attributes to ValueWithDOI
        scalar_params = {
            key: ValueWithDOI(value) for key, value in scalar_params.items()
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
    startdls: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0),
        description="Start of daylight savings time in decimal day of year",
    )
    enddls: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0),
        description="End of daylight savings time in decimal day of year",
    )
    heat: AnthropogenicHeat = Field(
        description="Anthropogenic heat emission parameters",
        default_factory=AnthropogenicHeat,
    )
    co2: CO2Params = Field(
        description="CO2 emission parameters", default_factory=CO2Params
    )

    ref: Optional[Reference] = None

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
        df_state.loc[grid_id, ("startdls", "0")] = self.startdls.value
        df_state.loc[grid_id, ("enddls", "0")] = self.enddls.value

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
        startdls = ValueWithDOI(df.loc[grid_id, ("startdls", "0")])
        enddls = ValueWithDOI(df.loc[grid_id, ("enddls", "0")])

        # Reconstruct heat parameters
        heat = AnthropogenicHeat.from_df_state(df, grid_id)

        # Reconstruct CO2 parameters
        co2 = CO2Params.from_df_state(df, grid_id)

        return cls(startdls=startdls, enddls=enddls, heat=heat.model_dump(), co2=co2.model_dump())

class AnthropogenicHeat(
    BaseModel
):  # TODO: May need to add the ValueWithDOI to the profiles here
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

    ref: Optional[Reference] = None

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


class CO2Params(
    BaseModel
):  # TODO: May need to add the ValueWithDOI to the profiles here
    co2pointsource: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="CO2 point source emission factor"
    )
    ef_umolco2perj: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="CO2 emission factor per unit of fuel"
    )
    enef_v_jkm: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0),
        description="CO2 emission factor per unit of vehicle distance",
    )
    fcef_v_kgkm: DayProfile = Field(
        description="Fuel consumption efficiency for vehicles",
        default_factory=DayProfile,
    )
    frfossilfuel_heat: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Fraction of fossil fuel heat"
    )
    frfossilfuel_nonheat: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Fraction of fossil fuel non-heat"
    )
    maxfcmetab: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Maximum fuel consumption metabolic rate"
    )
    maxqfmetab: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Maximum heat production metabolic rate"
    )
    minfcmetab: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Minimum fuel consumption metabolic rate"
    )
    minqfmetab: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Minimum heat production metabolic rate"
    )
    trafficrate: DayProfile = Field(
        description="Traffic rate", default_factory=DayProfile
    )
    trafficunits: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0), description="Traffic units"
    )
    traffprof_24hr: HourlyProfile = Field(
        description="24-hour profile of traffic rate", default_factory=HourlyProfile
    )
    humactivity_24hr: HourlyProfile = Field(
        description="24-hour profile of human activity", default_factory=HourlyProfile
    )

    ref: Optional[Reference] = None

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
            "co2pointsource": self.co2pointsource.value,
            "ef_umolco2perj": self.ef_umolco2perj.value,
            "enef_v_jkm": self.enef_v_jkm.value,
            "frfossilfuel_heat": self.frfossilfuel_heat.value,
            "frfossilfuel_nonheat": self.frfossilfuel_nonheat.value,
            "maxfcmetab": self.maxfcmetab.value,
            "maxqfmetab": self.maxqfmetab.value,
            "minfcmetab": self.minfcmetab.value,
            "minqfmetab": self.minqfmetab.value,
            "trafficunits": self.trafficunits.value,
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

        # Convert scalar attributes to ValueWithDOI
        scalar_params = {
            key: ValueWithDOI(value) for key, value in scalar_params.items()
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

