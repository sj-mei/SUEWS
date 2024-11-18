from typing import Dict, List, Optional, Union, Literal, Tuple
from pydantic import (
    BaseModel,
    Field,
    model_validator,
    field_validator,
    validator,
    PrivateAttr,
)
import numpy as np
from enum import Enum
import pandas as pd
import yaml
import pdb


class SurfaceType(str, Enum):
    PAVED = "paved"
    BLDGS = "bldgs"
    DECTR = "dectr"
    EVETR = "evetr"
    GRASS = "grass"
    BSOIL = "bsoil"
    WATER = "water"


class SnowAlb(BaseModel):
    snowalb: float = Field(ge=0, le=1, description="Snow albedo")


class SurfaceInitialState(BaseModel):
    """Initial state parameters for a surface type"""

    state: float = Field(ge=0, description="Initial state of the surface")
    soilstore: float = Field(ge=0, description="Initial soil store")
    snowfrac: Optional[float] = Field(None, ge=0, le=1, description="Snow fraction")
    snowpack: Optional[float] = Field(None, ge=0, description="Snow pack")
    snowwater: Optional[float] = Field(None, ge=0, description="Snow water")
    snowdens: Optional[float] = Field(None, ge=0, description="Snow density")
    alb_id: Optional[float] = Field(
        None, description="Initial albedo for vegetated surfaces"
    )
    porosity_id: Optional[float] = Field(
        None, description="Initial porosity for deciduous trees"
    )
    decidcap_id: Optional[float] = Field(
        None, description="Initial deciduous capacity for deciduous trees"
    )
    lai_id: Optional[float] = Field(
        None, description="Initial leaf area index for vegetated surfaces"
    )

    _surface_type: Optional[SurfaceType] = PrivateAttr(default=None)

    def set_surface_type(self, surface_type: SurfaceType):
        """Set and validate surface type specific parameters"""
        self._surface_type = surface_type

    @model_validator(mode="after")
    def validate_surface_state(self) -> "SurfaceInitialState":
        """Validate state based on surface type if set"""
        if self._surface_type is not None:
            # Re-run validation with current surface type
            self.set_surface_type(self._surface_type)
        # Validate vegetation-specific parameters
        veg_types = [SurfaceType.DECTR, SurfaceType.EVETR, SurfaceType.GRASS]
        if self._surface_type in veg_types:
            if self.alb_id is None:
                raise ValueError(f"alb_id is required for {self._surface_type.value}")
            if self.lai_id is None:
                raise ValueError(f"lai_id is required for {self._surface_type.value}")

        # Validate deciduous-specific parameters
        if self._surface_type == SurfaceType.DECTR:
            if self.decidcap_id is None:
                raise ValueError("decidcap_id is required for deciduous trees")
            if self.porosity_id is None:
                raise ValueError("porosity_id is required for deciduous trees")

        # Validate non-vegetation parameters
        non_veg_types = [
            SurfaceType.PAVED,
            SurfaceType.BLDGS,
            SurfaceType.BSOIL,
            SurfaceType.WATER,
        ]
        if self._surface_type in non_veg_types:
            if any(
                param is not None
                for param in [
                    self.alb_id,
                    self.lai_id,
                    self.decidcap_id,
                    self.porosity_id,
                ]
            ):
                raise ValueError(
                    f"Vegetation parameters should not be set for {self._surface_type.value}"
                )
        return self


class InitialStates(BaseModel):
    """Initial conditions for the SUEWS model"""

    snowalb: float = Field(ge=0, le=1, description="Initial snow albedo")
    paved: SurfaceInitialState
    bldgs: SurfaceInitialState
    dectr: SurfaceInitialState
    evetr: SurfaceInitialState
    grass: SurfaceInitialState
    bsoil: SurfaceInitialState
    water: SurfaceInitialState

    def __init__(self, **data):
        super().__init__(**data)
        # Set surface types for each surface
        self.paved.set_surface_type(SurfaceType.PAVED)
        self.bldgs.set_surface_type(SurfaceType.BLDGS)
        self.dectr.set_surface_type(SurfaceType.DECTR)
        self.evetr.set_surface_type(SurfaceType.EVETR)
        self.grass.set_surface_type(SurfaceType.GRASS)
        self.bsoil.set_surface_type(SurfaceType.BSOIL)
        self.water.set_surface_type(SurfaceType.WATER)


class ThermalLayer(BaseModel):
    dz: List[float] = Field(min_items=5, max_items=5)
    k: List[float] = Field(min_items=5, max_items=5)
    cp: List[float] = Field(min_items=5, max_items=5)
    temperature: List[float] = Field(min_items=5, max_items=5)


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


class StorageDrainParams(BaseModel):
    store_min: float = Field(ge=0)
    store_max: float = Field(ge=0)
    store_cap: float = Field(ge=0)
    drain_eq: int
    drain_coef_1: float
    drain_coef_2: float


class OHMCoefficients(BaseModel):
    a1: Dict[str, float]
    a2: Dict[str, float]
    a3: Dict[str, float]


class SurfaceProperties(BaseModel):
    """Base properties for all surface types"""

    sfr: float = Field(ge=0, le=1, description="Surface fraction")
    emis: float = Field(ge=0, le=1, description="Surface emissivity")
    chanohm: Optional[float] = None
    cpanohm: Optional[float] = None
    kkanohm: Optional[float] = None
    ohm_threshsw: Optional[float] = None
    ohm_threshwd: Optional[float] = None
    ohm_coef: Optional[OHMCoefficients] = None
    soildepth: float
    soilstorecap: float
    statelimit: float
    sathydraulicconduct: float
    waterdist: Optional[WaterDistribution] = None
    storedrainprm: StorageDrainParams
    snowpacklimit: float = Field(ge=0, description="Snow pack limit for the surface")
    _surface_type: Optional[SurfaceType] = PrivateAttr(default=None)

    def set_surface_type(self, surface_type: SurfaceType):
        self._surface_type = surface_type
        # Validate water distribution after setting surface type
        if self._surface_type == SurfaceType.WATER:
            if self.waterdist is not None:
                raise ValueError("Water surface should not have water distribution")
        else:
            if self.waterdist is None:
                raise ValueError(
                    f"Water distribution required for {self._surface_type.value}"
                )
            self.waterdist.validate_distribution(self._surface_type)


class NonVegetatedSurfaceProperties(SurfaceProperties):
    alb: float = Field(ge=0, le=1, description="Surface albedo")


class PavedProperties(NonVegetatedSurfaceProperties):
    surface_type: Literal[SurfaceType.PAVED] = SurfaceType.PAVED


class BuildingProperties(NonVegetatedSurfaceProperties):
    surface_type: Literal[SurfaceType.BLDGS] = SurfaceType.BLDGS
    faibldg: float = Field(ge=0, description="Frontal area index of buildings")
    bldgh: float = Field(ge=0, description="Building height")


class BaresoilProperties(NonVegetatedSurfaceProperties):
    surface_type: Literal[SurfaceType.BSOIL] = SurfaceType.BSOIL


class WaterProperties(NonVegetatedSurfaceProperties):
    surface_type: Literal[SurfaceType.WATER] = SurfaceType.WATER
    flowchange: float
    irrfracwater: float


class ModelControl(BaseModel):
    tstep: int = Field(description="Time step in seconds")
    forcing_file: str
    output_file: str
    daylightsaving_method: int


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
    snowuse: int


class LUMPSParams(BaseModel):
    raincover: float = Field(ge=0, le=1)
    rainmaxres: float = Field(ge=0, le=1)
    drainrt: float = Field(ge=0, le=1)
    veg_type: int


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


class DayProfile(BaseModel):
    working_day: float
    holiday: float


class WeeklyProfile(BaseModel):
    monday: float
    tuesday: float
    wednesday: float
    thursday: float
    friday: float
    saturday: float
    sunday: float


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


class CO2Params(BaseModel):
    co2pointsource: float
    ef_umolco2perj: float
    enef_v_jkm: float
    fcef_v_kgkm: DayProfile
    frfossilfuel_heat: float
    frfossilfuel_nonheat: float
    maxfcmetab: float
    maxqfmetab: float
    min_res_bioco2: float = Field(default=0.1, description="Minimum respiratory biogenic CO2")
    minfcmetab: float
    minqfmetab: float
    trafficrate: DayProfile
    trafficunits: float
    traffprof_24hr: HourlyProfile
    humactivity_24hr: HourlyProfile


class AnthropogenicEmissions(BaseModel):
    startdls: float
    enddls: float
    heat: AnthropogenicHeat
    co2: CO2Params


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


class LAIPowerCoefficients(BaseModel):
    growth_lai: float = Field(
        description="Power coefficient for LAI in growth equation (LAIPower[1])"
    )
    growth_gdd: float = Field(
        description="Power coefficient for GDD in growth equation (LAIPower[2])"
    )
    senescence_lai: float = Field(
        description="Power coefficient for LAI in senescence equation (LAIPower[3])"
    )
    senescence_sdd: float = Field(
        description="Power coefficient for SDD in senescence equation (LAIPower[4])"
    )

    def to_list(self) -> List[float]:
        """Convert to list format for Fortran interface"""
        return [
            self.growth_lai,
            self.growth_gdd,
            self.senescence_lai,
            self.senescence_sdd,
        ]


class LAIParams(BaseModel):
    baset: float = Field(
        description="Base Temperature for initiating growing degree days (GDD) for leaf growth [degC]"
    )
    gddfull: float = Field(
        description="Growing degree days (GDD) needed for full capacity of LAI [degC]"
    )
    basete: float = Field(
        description="Base temperature for initiating senescence degree days (SDD) for leaf off [degC]"
    )
    sddfull: float = Field(
        description="Senescence degree days (SDD) needed to initiate leaf off [degC]"
    )
    laimin: float = Field(description="Leaf-off wintertime value [m2 m-2]")
    laimax: float = Field(description="Full leaf-on summertime value [m2 m-2]")
    laipower: LAIPowerCoefficients = Field(
        description="LAI calculation power parameters for growth and senescence"
    )
    laitype: int = Field(
        description="LAI calculation choice (0: original, 1: new high latitude)"
    )

    @model_validator(mode="after")
    def validate_lai_ranges(self) -> "LAIParams":
        if self.laimin > self.laimax:
            raise ValueError("laimin must be less than or equal to laimax")
        if self.baset > self.gddfull:
            raise ValueError("baset must be less than gddfull")
        return self


class VegetatedSurfaceProperties(SurfaceProperties):
    alb_min: float = Field(ge=0, le=1, description="Minimum albedo")
    alb_max: float = Field(ge=0, le=1, description="Maximum albedo")
    beta_bioco2: float
    beta_enh_bioco2: float
    alpha_bioco2: float
    alpha_enh_bioco2: float
    resp_a: float
    resp_b: float
    theta_bioco2: float
    maxconductance: float = Field(default=0.5, description="Maximum surface conductance")
    min_res_bioco2: float = Field(default=0.1, description="Minimum respiratory biogenic CO2")
    lai: LAIParams
    ie_a: float = Field(description="Irrigation efficiency coefficient-automatic")
    ie_m: float = Field(description="Irrigation efficiency coefficient-manual")

    @model_validator(mode="after")
    def validate_albedo_range(self) -> "VegetatedSurfaceProperties":
        if self.alb_min > self.alb_max:
            error_message = ValueError(
                f"alb_min (input {self.alb_min}) must be less than or equal to alb_max (entered {self.alb_max})."
            )
            exceptions.append(error_message)
            # raise ValueError(f"alb_min (input {self.alb_min}) must be less than or equal to alb_max (entered {self.alb_max}).")
        return self


class DectrProperties(VegetatedSurfaceProperties):
    faidectree: float = Field(description="Frontal area index of deciduous trees")
    dectreeh: float = Field(description="Height of deciduous trees [m]")
    pormin_dec: float = Field(
        ge=0.1,
        le=0.9,
        description="Minimum porosity of deciduous trees in winter when leaves are off",
    )
    pormax_dec: float = Field(
        ge=0.1,
        le=0.9,
        description="Maximum porosity of deciduous trees in summer when leaves are fully on",
    )
    capmax_dec: float = Field(
        description="Maximum storage capacity of deciduous trees in summer when leaves are fully on [mm]"
    )
    capmin_dec: float = Field(
        description="Minimum storage capacity of deciduous trees in winter when leaves are off [mm]"
    )

    @model_validator(mode="after")
    def validate_porosity_range(self) -> "DectrProperties":
        if self.pormin_dec >= self.pormax_dec:
            error_message = ValueError(
                f"pormin_dec ({self.pormin_dec}) must be less than pormax_dec ({self.pormax_dec})."
            )
            exceptions.append(error_message)
            # raise ValueError(f"pormin_dec ({self.pormin_dec}) must be less than pormax_dec ({self.pormax_dec}).")
        return self

    @model_validator(mode="after")
    def validate_cap_range(self) -> "DectrProperties":
        if self.capmin_dec >= self.capmax_dec:
            error_message = ValueError(
                f"capmin_dec ({self.capmin_dec}) must be less than capmax_dec ({self.capmax_dec})."
            )
            exceptions.append(error_message)
            # raise ValueError(f"capmin_dec ({self.capmin_dec}) must be less than capmax_dec ({self.capmax_dec}).")
        return self


class EvetrProperties(VegetatedSurfaceProperties):
    faievetree: float = Field(description="Frontal area index of evergreen trees")
    evetreeh: float = Field(description="Height of evergreen trees [m]")


class SnowParams(BaseModel):
    crwmax: float
    crwmin: float
    narp_emis_snow: float
    preciplimit: float
    preciplimitalb: float
    snowalbmax: float
    snowalbmin: float
    snowdensmin: float
    snowlimbldg: float
    snowlimpaved: float
    snowprof_24hr: HourlyProfile
    tau_a: float
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


class Site(BaseModel):
    name: str
    properties: SiteProperties
    initial_states: InitialStates

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
        df = pd.DataFrame(index=[0], columns=columns)

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
            df.loc[0, (col_name, str_indices)] = value

        # Model control
        set_df_value("tstep", 0, self.model.control.tstep)

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
        set_df_value("snowuse", 0, self.model.physics.snowuse)


        # Process each site (assuming single site for now)
        site = self.site[0]
        props = site.properties

        # Basic site properties
        set_df_value("lat", 0, props.lat)
        set_df_value("lng", 0, props.lng)
        set_df_value("alt", 0, props.alt)
        set_df_value("timezone", 0, props.timezone)
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
        set_df_value("ground_albedo_dir_mult_fact", 0, spartacus.ground_albedo_dir_mult_fact)
        set_df_value("n_stream_lw_urban", 0, spartacus.n_stream_lw_urban)
        set_df_value("n_stream_sw_urban", 0, spartacus.n_stream_sw_urban)
        set_df_value("n_vegetation_region_urban", 0, spartacus.n_vegetation_region_urban)
        set_df_value("sw_dn_direct_frac", 0, spartacus.sw_dn_direct_frac)
        set_df_value("use_sw_direct_albedo", 0, spartacus.use_sw_direct_albedo)
        set_df_value("veg_contact_fraction_const", 0, spartacus.veg_contact_fraction_const)
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
            set_df_value("wuprofa_24hr", (hour, 0), irr.wuprofa_24hr.working_day[hour_str])
            set_df_value("wuprofa_24hr", (hour, 1), irr.wuprofa_24hr.holiday[hour_str])

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

            # Basic surface properties
            set_df_value("sfr_surf", (surf_idx,), surface.sfr)
            set_df_value("emis", (surf_idx,), surface.emis)

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

            # Handle albedo
            if hasattr(surface, "alb"):
                set_df_value("alb", (surf_idx,), surface.alb)
            elif hasattr(surface, "alb_min") and hasattr(surface, "alb_max"):
                set_df_value("alb", (surf_idx,), surface.alb_min)  # Use min as default
                if surf_name == "dectr":
                    set_df_value("albmax_dectr", 0, surface.alb_max)
                    set_df_value("albmin_dectr", 0, surface.alb_min)
                elif surf_name == "evetr":
                    set_df_value("albmax_evetr", 0, surface.alb_max)
                    set_df_value("albmin_evetr", 0, surface.alb_min)
                elif surf_name == "grass":
                    set_df_value("albmax_grass", 0, surface.alb_max)
                    set_df_value("albmin_grass", 0, surface.alb_min)
            # Add to OHM section
            if hasattr(surface, "ohm_threshsw"):
                set_df_value("ohm_threshsw", (surf_idx,), surface.ohm_threshsw)
            if hasattr(surface, "ohm_threshwd"):
                set_df_value("ohm_threshwd", (surf_idx,), surface.ohm_threshwd)
            if hasattr(surface, "chanohm"):
                set_df_value("chanohm", (surf_idx,), surface.chanohm)
            if hasattr(surface, "cpanohm"):
                set_df_value("cpanohm", (surf_idx,), surface.cpanohm)
            if hasattr(surface, "kkanohm"):
                set_df_value("kkanohm", (surf_idx,), surface.kkanohm)

            # Add to surface properties section
            if hasattr(surface, "soildepth"):
                set_df_value("soildepth", (surf_idx,), surface.soildepth)
            if hasattr(surface, "soilstorecap"):
                set_df_value("soilstorecap_surf", (surf_idx,), surface.soilstorecap)
            if hasattr(surface, "statelimit"):
                set_df_value("statelimit_surf", (surf_idx,), surface.statelimit)

            # Initial states
            init_state = getattr(site.initial_states, surf_name)
            if init_state:
                # Basic state parameters
                set_df_value("state_surf", (surf_idx,), init_state.state)
                set_df_value("soilstore_surf", (surf_idx,), init_state.soilstore)

                # Snow-related parameters
                if init_state.snowfrac is not None:
                    set_df_value("snowfrac", (surf_idx,), init_state.snowfrac)
                if init_state.snowpack is not None:
                    set_df_value("snowpack", (surf_idx,), init_state.snowpack)
                if init_state.snowwater is not None:
                    set_df_value("snowwater", (surf_idx,), init_state.snowwater)
                if init_state.snowdens is not None:
                    set_df_value("snowdens", (surf_idx,), init_state.snowdens)

                # Vegetation-specific parameters
                if surf_name in ["dectr", "evetr", "grass"]:
                    # albedo
                    if init_state.alb_id is not None:
                        set_df_value(f"alb{surf_name}_id", 0, init_state.alb_id)
                    # LAI
                    if init_state.lai_id is not None:
                        set_df_value(f"lai_id", (surf_idx - 2,), init_state.lai_id)
                    # Additional parameters for deciduous trees
                    if surf_name == "dectr":
                        if init_state.decidcap_id is not None:
                            set_df_value("decidcap_id", 0, init_state.decidcap_id)
                        if init_state.porosity_id is not None:
                            set_df_value("porosity_id", 0, init_state.porosity_id)

                    # Missing porosity parameters for deciduous trees
                    if hasattr(surface, "pormin_dec"):
                        set_df_value("pormin_dec", 0, surface.pormin_dec)
                    if hasattr(surface, "pormax_dec"):
                        set_df_value("pormax_dec", 0, surface.pormax_dec)

                    # Missing capacity parameters for deciduous trees
                    if hasattr(surface, "capmin_dec"):
                        set_df_value("capmin_dec", 0, surface.capmin_dec)
                    if hasattr(surface, "capmax_dec"):
                        set_df_value("capmax_dec", 0, surface.capmax_dec)
                    if hasattr(surface, "min_res_bioco2"):
                        set_df_value("min_res_bioco2", (surf_idx - 2,), surface.min_res_bioco2)

            # Set initial snow albedo
            set_df_value("snowalb", 0, site.initial_states.snowalb)

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
                prop_day = getattr(irrigation.wuprofa_24hr, day)
                set_df_value(f"wuprofa_24hr", (hour, i), prop_day[f"{hour+1}"])
                prop_day = getattr(irrigation.wuprofm_24hr, day)
                set_df_value(f"wuprofm_24hr", (hour, i), prop_day[f"{hour+1}"])

        # Process anthropogenic emissions
        anthro = props.anthropogenic_emissions
        set_df_value("startdls", 0, anthro.startdls)
        set_df_value("enddls", 0, anthro.enddls)
        # heat
        for i, day in enumerate(["working_day", "holiday"]):
            set_df_value("qf0_beu", (i,), getattr(anthro.heat.qf0_beu, day))
            set_df_value("qf_a", (i,), getattr(anthro.heat.qf_a, day))
            set_df_value("qf_b", (i,), getattr(anthro.heat.qf_b, day))
            set_df_value("qf_c", (i,), getattr(anthro.heat.qf_c, day))
            set_df_value("baset_cooling", (i,), getattr(anthro.heat.baset_cooling, day))
            set_df_value("baset_heating", (i,), getattr(anthro.heat.baset_heating, day))
            set_df_value("ah_min", (i,), getattr(anthro.heat.ah_min, day))
            set_df_value(
                "ah_slope_cooling", (i,), getattr(anthro.heat.ah_slope_cooling, day)
            )
            set_df_value(
                "ah_slope_heating", (i,), getattr(anthro.heat.ah_slope_heating, day)
            )
            # Hourly profiles
            for hour in range(24):
                set_df_value(
                    "ahprof_24hr",
                    (hour, i),
                    getattr(anthro.heat.ahprof_24hr, day)[f"{hour+1}"],
                )

        # CO2 parameters
        set_df_value("co2pointsource", 0, anthro.co2.co2pointsource)
        set_df_value("ef_umolco2perj", 0, anthro.co2.ef_umolco2perj)
        set_df_value("enef_v_jkm", 0, anthro.co2.enef_v_jkm)

        # Traffic and emission factors for weekday/weekend
        for i, day in enumerate(["working_day", "holiday"]):
            # Traffic rate and emission factors
            set_df_value("trafficrate", (i,), getattr(anthro.co2.trafficrate, day))
            set_df_value("fcef_v_kgkm", (i,), getattr(anthro.co2.fcef_v_kgkm, day))

            # 24-hour profiles
            for hour in range(24):
                set_df_value(
                    "traffprof_24hr",
                    (hour, i),
                    getattr(anthro.co2.traffprof_24hr, day)[f"{hour+1}"],
                )
                set_df_value(
                    "humactivity_24hr",
                    (hour, i),
                    getattr(anthro.co2.humactivity_24hr, day)[f"{hour+1}"],
                )

        # Traffic units setting
        set_df_value("trafficunits", 0, anthro.co2.trafficunits)

        # Additional emission-related parameters
        set_df_value("frfossilfuel_heat", 0, anthro.co2.frfossilfuel_heat)
        set_df_value("frfossilfuel_nonheat", 0, anthro.co2.frfossilfuel_nonheat)
        set_df_value("maxfcmetab", 0, anthro.co2.maxfcmetab)
        set_df_value("maxqfmetab", 0, anthro.co2.maxqfmetab)
        set_df_value("minfcmetab", 0, anthro.co2.minfcmetab)
        set_df_value("minqfmetab", 0, anthro.co2.minqfmetab)

        # Snow parameters
        snow = props.snow
        set_df_value("crwmax", 0, snow.crwmax)
        set_df_value("crwmin", 0, snow.crwmin)
        set_df_value("narp_emis_snow", 0, snow.narp_emis_snow)
        set_df_value("preciplimit", 0, snow.preciplimit)
        set_df_value("preciplimitalb", 0, snow.preciplimitalb)
        set_df_value("snowalbmax", 0, snow.snowalbmax)
        set_df_value("snowalbmin", 0, snow.snowalbmin)
        set_df_value("snowdensmin", 0, snow.snowdensmin)
        set_df_value("snowlimbldg", 0, snow.snowlimbldg)
        set_df_value("snowlimpaved", 0, snow.snowlimpaved)
        set_df_value("tau_a", 0, snow.tau_a)
        set_df_value("tau_r", 0, snow.tau_r)
        set_df_value("tempmeltfact", 0, snow.tempmeltfact)
        set_df_value("radmeltfact", 0, snow.radmeltfact)

        # Missing height parameters for vegetation and buildings
        set_df_value("bldgh", 0, props.land_cover.bldgs.bldgh)  # Building height
        set_df_value(
            "evetreeh", 0, props.land_cover.evetr.evetreeh
        )  # Evergreen tree height
        set_df_value(
            "dectreeh", 0, props.land_cover.dectr.dectreeh
        )  # Deciduous tree height
        # Add to to_df_state where vegetation parameters are handled


        # Missing FAI parameters
        set_df_value("faibldg", 0, props.land_cover.bldgs.faibldg)
        set_df_value("faievetree", 0, props.land_cover.evetr.faievetree)
        set_df_value("faidectree", 0, props.land_cover.dectr.faidectree)

        # Add to surface properties section
        if hasattr(surface, "sathydraulicconduct"):
            set_df_value(
                "sathydraulicconduct", (surf_idx,), surface.sathydraulicconduct
            )

        # Irrigation fraction parameters
        if hasattr(surface, "irrfracevetr"):
            set_df_value("irrfracevetr", 0, surface.irrfracevetr)
        if hasattr(surface, "irrfracbsoil"):
            set_df_value("irrfracbsoil", 0, surface.irrfracbsoil)
        if hasattr(surface, "irrfracwater"):
            set_df_value("irrfracwater", 0, surface.irrfracwater)

        # Water specific parameters
        if hasattr(surface, "flowchange"):
            set_df_value("flowchange", 0, surface.flowchange)

        # LUMPS parameters
        lumps = props.lumps
        set_df_value("raincover", 0, lumps.raincover)
        set_df_value("rainmaxres", 0, lumps.rainmaxres)
        set_df_value("drainrt", 0, lumps.drainrt)
        set_df_value("veg_type", 0, lumps.veg_type)

        # Population profile parameters
        set_df_value("popdensnighttime", 0, anthro.heat.popdensnighttime)
        for i, day in enumerate(["working_day", "holiday"]):
            set_df_value("popdensdaytime", (i,), getattr(anthro.heat.popdensdaytime, day))
            # 24-hour population profile
            for hour in range(24):
                set_df_value(
                    "popprof_24hr",
                    (hour, i),
                    getattr(anthro.heat.popprof_24hr, day)[f"{hour+1}"],
                )

        # Additional CO2 parameters
        set_df_value("frfossilfuel_heat", 0, anthro.co2.frfossilfuel_heat)
        set_df_value("frfossilfuel_nonheat", 0, anthro.co2.frfossilfuel_nonheat)
        set_df_value("maxfcmetab", 0, anthro.co2.maxfcmetab)
        set_df_value("maxqfmetab", 0, anthro.co2.maxqfmetab)
        set_df_value("minfcmetab", 0, anthro.co2.minfcmetab)
        set_df_value("minqfmetab", 0, anthro.co2.minqfmetab)
        set_df_value("trafficrate", (0,), anthro.co2.trafficrate.working_day)
        set_df_value("trafficrate", (1,), anthro.co2.trafficrate.holiday)
        set_df_value("trafficunits", 0, anthro.co2.trafficunits)

        # Snow profile
        for i, day in enumerate(["working_day", "holiday"]):
            for hour in range(24):
                set_df_value(
                    "snowprof_24hr",
                    (hour, i),
                    getattr(snow.snowprof_24hr, day)[f"{hour+1}"],
                )

        # Add to CO2 parameters section for vegetated surfaces
        if hasattr(surface, "resp_a"):
            set_df_value("resp_a", (surf_idx - 2,), surface.resp_a)
        if hasattr(surface, "resp_b"):
            set_df_value("resp_b", (surf_idx - 2,), surface.resp_b)
        if hasattr(surface, "theta_bioco2"):
            set_df_value("theta_bioco2", (surf_idx - 2,), surface.theta_bioco2)

        return df

    @classmethod
    def from_df_state(cls, df: pd.DataFrame) -> "SUEWSConfig":
        """Create config from DataFrame state"""
        surface_map = {
            0: "paved",
            1: "bldgs",
            2: "dectr",
            3: "evetr",
            4: "grass",
            5: "bsoil",
            6: "water",
        }

        # Initialize basic structure
        config_dict = {
            "name": "Generated from df_state",
            "description": "Automatically converted from DataFrame state",
            "model": {
                "control": {
                    "tstep": 300,  # Default values
                    "forcing_file": "forcing.csv",
                    "output_file": "output.csv",
                    "daylightsaving_method": 1,
                },
                "physics": {
                    # Add physics parameters from df
                },
            },
            "site": [
                {
                    "name": "site_0",
                    "properties": {
                        "lat": float(df.loc[0, ("lat", 0)]),
                        "lng": float(df.loc[0, ("lng", 0)]),
                        "alt": float(df.loc[0, ("alt", 0)]),
                        "timezone": int(df.loc[0, ("timezone", 0)]),
                        # ... other properties ...
                    },
                }
            ],
        }

        # Convert DataFrame data to config structure
        for surf_idx, surf_name in surface_map.items():
            # Extract surface properties
            surface_props = {
                "sfr": float(df.loc[0, ("sfr_surf", surf_idx)]),
                "emis": float(df.loc[0, ("emis", surf_idx)]),
                # ... other properties ...
            }

            # Add to config dictionary
            config_dict["site"][0]["properties"]["land_cover"][surf_name] = (
                surface_props
            )

        return cls(**config_dict)


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


# # Convert back to config
# suews_config_back = SUEWSConfig.from_df_state(df_state)
# print("testing from_df_state done!")
