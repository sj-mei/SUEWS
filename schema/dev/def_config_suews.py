from typing import Dict, List, Optional, Union, Literal, Tuple
from pydantic import BaseModel, Field, model_validator, field_validator, validator, PrivateAttr
import numpy as np
from enum import Enum
import pandas as pd
import yaml


class SurfaceType(str, Enum):
    PAVED = "paved"
    BLDGS = "bldgs"
    DECTR = "dectr"
    EVETR = "evetr"
    GRASS = "grass"
    BSOIL = "bsoil"
    WATER = "water"


class InitialState(BaseModel):
    state: float = Field(ge=0, description="Initial state")
    soilstore: float = Field(ge=0, description="Initial soil store")
    snowfrac: float = Field(ge=0, le=1, description="Snow fraction")
    snowpack: float = Field(ge=0, description="Snow pack")
    snowwater: float = Field(ge=0, description="Snow water")
    snowalb: float = Field(ge=0, le=1, description="Snow albedo")
    snowdens: float = Field(ge=0, description="Snow density")

    # Additional fields for specific surface types
    decidcap_id: Optional[float] = Field(None, description="Deciduous capacity ID (required for deciduous trees)")
    porosity_id: Optional[float] = Field(None, description="Porosity ID (required for deciduous trees)")
    alb_id: Optional[float] = Field(None, description="Albedo ID (required for vegetated surfaces)")
    surface_type: Optional[SurfaceType] = Field(None, description="Surface type for validation")

    @model_validator(mode='after')
    def validate_surface_specific_fields(self) -> 'InitialState':
        if self.surface_type == SurfaceType.DECTR:
            if self.decidcap_id is None:
                raise ValueError("decidcap_id is required for deciduous trees")
            if self.porosity_id is None:
                raise ValueError("porosity_id is required for deciduous trees")

        veg_types = [SurfaceType.DECTR, SurfaceType.EVETR, SurfaceType.GRASS]
        if self.surface_type in veg_types and self.alb_id is None:
            raise ValueError("alb_id is required for vegetated surfaces")

        return self


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
            SurfaceType.PAVED: ["to_bldgs", "to_dectr", "to_evetr", "to_grass", "to_bsoil", "to_water", "to_runoff"],
            SurfaceType.BLDGS: ["to_paved", "to_dectr", "to_evetr", "to_grass", "to_bsoil", "to_water", "to_runoff"],
            SurfaceType.DECTR: ["to_paved", "to_bldgs", "to_evetr", "to_grass", "to_bsoil", "to_water", "to_soilstore"],
            SurfaceType.EVETR: ["to_paved", "to_bldgs", "to_dectr", "to_grass", "to_bsoil", "to_water", "to_soilstore"],
            SurfaceType.GRASS: ["to_paved", "to_bldgs", "to_dectr", "to_evetr", "to_bsoil", "to_water", "to_soilstore"],
            SurfaceType.BSOIL: ["to_paved", "to_bldgs", "to_dectr", "to_evetr", "to_grass", "to_water", "to_soilstore"],
            SurfaceType.WATER: None  # Water surface doesn't have water distribution
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
    _surface_type: Optional[SurfaceType] = PrivateAttr(default=None)  # Private attribute for surface type

    def set_surface_type(self, surface_type: SurfaceType):
        self._surface_type = surface_type
        # Validate water distribution after setting surface type
        if self._surface_type == SurfaceType.WATER:
            if self.waterdist is not None:
                raise ValueError("Water surface should not have water distribution")
        else:
            if self.waterdist is None:
                raise ValueError(f"Water distribution required for {self._surface_type.value}")
            self.waterdist.validate_distribution(self._surface_type)


class NonVegetatedSurfaceProperties(SurfaceProperties):
    alb: float = Field(ge=0, le=1, description="Surface albedo")


class PavedProperties(NonVegetatedSurfaceProperties):
    surface_type: Literal[SurfaceType.PAVED] = SurfaceType.PAVED


class BuildingProperties(NonVegetatedSurfaceProperties):
    surface_type: Literal[SurfaceType.BLDGS] = SurfaceType.BLDGS
    faibldg: float = Field(ge=0, description="Frontal area index of buildings")


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
    aerodynamicresistancemethod: int
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

    @field_validator('working_day', 'holiday', mode='before')
    def convert_keys_to_str(cls, v: Dict) -> Dict[str, float]:
        if isinstance(v, dict):
            return {str(k): float(v) for k, v in v.items()}
        return v

    @model_validator(mode='after')
    def validate_hours(self) -> 'HourlyProfile':
        for profile in [self.working_day, self.holiday]:
            hours = [int(h) for h in profile.keys()]
            if not all(1 <= h <= 24 for h in hours):
                raise ValueError("Hour values must be between 1 and 24")
            if sorted(hours) != list(range(1, 25)):
                raise ValueError("Must have all hours from 1 to 24")
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
    ah_min: float
    ah_slope_cooling: float
    ah_slope_heating: float
    ahprof_24hr: HourlyProfile


class CO2Params(BaseModel):
    co2pointsource: float
    ef_umolco2perj: float
    enef_v_jkm: float
    fcef_v_kgkm: float
    frfossilfuel_heat: float
    frfossilfuel_nonheat: float
    maxfcmetab: float
    maxqfmetab: float
    min_res_bioco2: float
    minfcmetab: float
    minqfmetab: float
    trafficrate: float
    trafficunits: float
    traffprof_24hr: HourlyProfile


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
    maxconductance: float = Field(description="Maximum surface conductance")
    s1: float = Field(description="Soil moisture threshold parameter")
    s2: float = Field(description="Soil moisture threshold parameter")


class LAIParams(BaseModel):
    baset: float
    gddfull: float
    basete: float
    sddfull: float
    laimin: float
    laimax: float
    laipower: float
    laitype: int


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
    conductance: Conductance
    lai: LAIParams
    ie_a: float = Field(description="Irrigation efficiency coefficient-automatic")
    ie_m: float = Field(description="Irrigation efficiency coefficient-manual")

    @model_validator(mode='after')
    def validate_albedo_range(self) -> 'VegetatedSurfaceProperties':
        if self.alb_min > self.alb_max:
            raise ValueError("alb_min must be less than or equal to alb_max")
        return self


class DectrProperties(VegetatedSurfaceProperties):
    faidectree: float
    dectreeh: float
    pormin_dec: float
    pormax_dec: float
    capmax_dec: float
    capmin_dec: float


class EvetrProperties(VegetatedSurfaceProperties):
    faievetree: float
    evetreeh: float


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
    snowpacklimit: float
    snowprof_24hr: HourlyProfile
    tau_a: float
    tau_r: float
    tempmeltfact: float
    radmeltfact: float


class LandCover(BaseModel):
    paved: PavedProperties
    bldgs: BuildingProperties
    dectr: DectrProperties
    evetr: EvetrProperties
    grass: VegetatedSurfaceProperties
    bsoil: BaresoilProperties
    water: WaterProperties

    @model_validator(mode='after')
    def set_surface_types(self) -> 'LandCover':
        # Set surface types and validate
        surface_map = {
            'paved': (self.paved, SurfaceType.PAVED),
            'bldgs': (self.bldgs, SurfaceType.BLDGS),
            'dectr': (self.dectr, SurfaceType.DECTR),
            'evetr': (self.evetr, SurfaceType.EVETR),
            'grass': (self.grass, SurfaceType.GRASS),
            'bsoil': (self.bsoil, SurfaceType.BSOIL),
            'water': (self.water, SurfaceType.WATER),
        }

        for prop, surface_type in surface_map.values():
            prop.set_surface_type(surface_type)

        return self


class InitialConditions(BaseModel):
    soilstore_id: float
    wetstore_id: float
    snowstore_id: float
    snowwater_id: float
    snowpack_id: float
    snowdens_id: float
    snowfrac_id: float
    snowmelt_id: float
    snowalb_id: float
    gdd_id: float
    sdd_id: float
    laistore_id: float
    surface_temp_id: float
    surface_wetness_id: float
    surface_state_id: float
    surface_gc_id: float


class SiteProperties(BaseModel):
    lat: float = Field(ge=-90, le=90)
    lon: float = Field(ge=-180, le=180)
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
    irrigation: IrrigationParams
    anthropogenic_emissions: AnthropogenicEmissions
    snow: SnowParams
    land_cover: LandCover


class Site(BaseModel):
    name: str
    properties: SiteProperties
    initial_states: Dict[str, InitialState]


class SUEWSConfig(BaseModel):
    name: str
    description: str
    model: Dict[str, Union[ModelControl, ModelPhysics]]
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
            if indices == "0":
                tuples.append((col_name, 0))
            else:
                # Convert string tuple to actual tuple
                indices = eval(indices)
                if isinstance(indices, int):
                    indices = (indices,)
                tuples.append((col_name, *indices))

        return pd.MultiIndex.from_tuples(tuples)

    def to_df_state(self) -> pd.DataFrame:
        """Convert config to DataFrame state format"""
        # Initialize empty DataFrame with correct structure
        columns = self.create_multi_index_columns("df_state_columns.txt")
        df = pd.DataFrame(index=[0], columns=columns)

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, indices: Union[int, Tuple], value: float):
            if isinstance(indices, int):
                indices = (indices,)
            df.loc[0, (col_name, *indices)] = value

        # Process each site (assuming single site for now)
        site = self.site[0]
        props = site.properties

        # Basic site properties
        set_df_value('lat', 0, props.lat)
        set_df_value('lon', 0, props.lon)
        set_df_value('alt', 0, props.alt)
        set_df_value('timezone', 0, props.timezone)

        # Surface properties
        surface_map = {
            'paved': 0, 'bldgs': 1, 'dectr': 2, 'evetr': 3,
            'grass': 4, 'bsoil': 5, 'water': 6
        }

        # Process each surface type
        for surf_name, surf_idx in surface_map.items():
            surface = getattr(props.land_cover, surf_name)

            # Basic surface properties
            set_df_value('sfr', surf_idx, surface.sfr)
            set_df_value('alb', surf_idx,
                        getattr(surface, 'alb', (surface.alb_min + surface.alb_max) / 2))
            set_df_value('emis', surf_idx, surface.emis)

            # Initial states
            init_state = getattr(props.initial_states, surf_name)
            set_df_value('state_surf', surf_idx, init_state.state)
            set_df_value('soilstore_surf', surf_idx, init_state.soilstore)
            set_df_value('snowwater', surf_idx, init_state.snowwater)
            set_df_value('snowalb', surf_idx, init_state.snowalb)
            set_df_value('snowdens', surf_idx, init_state.snowdens)
            set_df_value('snowfrac', surf_idx, init_state.snowfrac)

            # Vegetation-specific properties
            if surf_name in ['dectr', 'evetr', 'grass']:
                if hasattr(init_state, 'alb_id'):
                    set_df_value(f'alb{surf_name}_id', 0, init_state.alb_id)
                if hasattr(init_state, 'porosity_id'):
                    set_df_value('porosity_id', 0, init_state.porosity_id)
                if hasattr(init_state, 'decidcap_id'):
                    set_df_value('decidcap_id', 0, init_state.decidcap_id)

        # Process anthropogenic emissions
        anthro = props.anthropogenic_emissions
        for hour in range(24):
            set_df_value('ahprof_24hr', (hour, 0),
                        anthro.heat.ahprof_24hr.working_day[str(hour + 1)])
            set_df_value('ahprof_24hr', (hour, 1),
                        anthro.heat.ahprof_24hr.holiday[str(hour + 1)])

        # ... Add more conversions as needed ...

        return df

    @classmethod
    def from_df_state(cls, df: pd.DataFrame) -> "SUEWSConfig":
        """Create config from DataFrame state"""
        surface_map = {
            0: 'paved', 1: 'bldgs', 2: 'dectr', 3: 'evetr',
            4: 'grass', 5: 'bsoil', 6: 'water'
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
                    "daylightsaving_method": 1
                },
                "physics": {
                    # Add physics parameters from df
                }
            },
            "site": [{
                "name": "site_0",
                "properties": {
                    "lat": float(df.loc[0, ('lat', 0)]),
                    "lon": float(df.loc[0, ('lon', 0)]),
                    "alt": float(df.loc[0, ('alt', 0)]),
                    "timezone": int(df.loc[0, ('timezone', 0)]),
                    # ... other properties ...
                }
            }]
        }

        # Convert DataFrame data to config structure
        for surf_idx, surf_name in surface_map.items():
            # Extract surface properties
            surface_props = {
                "sfr": float(df.loc[0, ('sfr', surf_idx)]),
                "emis": float(df.loc[0, ('emis', surf_idx)]),
                # ... other properties ...
            }

            # Add to config dictionary
            config_dict["site"][0]["properties"]["land_cover"][surf_name] = surface_props

        return cls(**config_dict)

# Create a sample configuration
# sample_config_suews = SUEWSConfig(
#     name="Sample Config",
#     description="A sample configuration for testing",
#     model={
#         "control": ModelControl(
#             tstep=300,
#             forcing_file="forcing.csv",
#             output_file="output.csv",
#             daylightsaving_method=1,
#         ),
#         "physics": ModelPhysics(
#             netradiationmethod=1,
#             emissionsmethod=1,
#             storageheatmethod=1,
#             ohmincqf=1,
#             aerodynamicresistancemethod=1,
#             roughlenmommethod=1,
#             roughlenheatmethod=1,
#             stabilitymethod=1,
#             smdmethod=1,
#             waterusemethod=1,
#             diagmethod=1,
#             snowuse=0,
#         ),
#     },
#     site=[
#         Site(
#             name="site_0",
#             properties=SiteProperties(
#                 lat=51.5,
#                 lon=-0.1,
#                 alt=35,
#                 timezone=0,
#                 surfacearea=1,
#                 z=10,
#                 z0m_in=0.1,
#                 zdm_in=2,
#                 pipecapacity=100,
#                 runofftowater=0.5,
#                 narp_trans_site=1,
#                 lumps=LUMPSParams(
#                     raincover=0.25,
#                     rainmaxres=0.25,
#                     drainrt=0.25,
#                     veg_type=1,
#                 ),
#                 spartacus=SPARTACUSParams(
#                     air_ext_lw=1.0,
#                     air_ext_sw=1.0,
#                     air_ssa_lw=1.0,
#                     air_ssa_sw=1.0,
#                     ground_albedo_dir_mult_fact=1.0,
#                     n_stream_lw_urban=1,
#                     n_stream_sw_urban=1,
#                     n_vegetation_region_urban=1,
#                     sw_dn_direct_frac=1.0,
#                     use_sw_direct_albedo=0.1,
#                     veg_contact_fraction_const=1.0,
#                     veg_fsd_const=1.0,
#                     veg_ssa_lw=1.0,
#                     veg_ssa_sw=1.0,
#                 ),
#                 irrigation=IrrigationParams(
#                     h_maintain=0.5,
#                     faut=0.5,
#                     ie_start=100,
#                     ie_end=200,
#                     internalwateruse_h=10,
#                     daywatper=DayProfile(
#                         Monday=10.0,
#                         Tuesday=10.0,
#                         Wednesday=10.0,
#                         Thursday=10.0,
#                         Friday=10.0,
#                         Saturday=10.0,
#                         Sunday=10.0,
#                     ),
#                     daywat=DayProfile(
#                         Monday=10.0,
#                         Tuesday=10.0,
#                         Wednesday=10.0,
#                         Thursday=10.0,
#                         Friday=10.0,
#                         Saturday=10.0,
#                         Sunday=10.0,
#                     ),
#                     wuprofa_24hr=HourlyProfile(
#                         working_day={str(i): 0.1 for i in range(1, 25)},
#                         holiday={str(i): 0.1 for i in range(1, 25)},
#                     ),
#                     wuprofm_24hr=HourlyProfile(
#                         working_day={str(i): 0.1 for i in range(1, 25)},
#                         holiday={str(i): 0.1 for i in range(1, 25)},
#                     ),
#                 ),
#                 anthropogenic_emissions=AnthropogenicEmissions(
#                     startdls=100,
#                     enddls=300,
#                     heat=AnthropogenicHeat(
#                         qf0_beu={
#                             "paved": 0.1,
#                             "bldgs": 0.1,
#                             "dectr": 0.1,
#                             "evetr": 0.1,
#                             "grass": 0.1,
#                             "bsoil": 0.1,
#                             "water": 0.1,
#                         },
#                         qf_a={
#                             "paved": 0.1,
#                             "bldgs": 0.1,
#                             "dectr": 0.1,
#                             "evetr": 0.1,
#                             "grass": 0.1,
#                             "bsoil": 0.1,
#                             "water": 0.1,
#                         },
#                         qf_b={
#                             "paved": 0.1,
#                             "bldgs": 0.1,
#                             "dectr": 0.1,
#                             "evetr": 0.1,
#                             "grass": 0.1,
#                             "bsoil": 0.1,
#                             "water": 0.1,
#                         },
#                         qf_c={
#                             "paved": 0.1,
#                             "bldgs": 0.1,
#                             "dectr": 0.1,
#                             "evetr": 0.1,
#                             "grass": 0.1,
#                             "bsoil": 0.1,
#                             "water": 0.1,
#                         },
#                         baset_cooling={
#                             "paved": 0.1,
#                             "bldgs": 0.1,
#                             "dectr": 0.1,
#                             "evetr": 0.1,
#                             "grass": 0.1,
#                             "bsoil": 0.1,
#                             "water": 0.1,
#                         },
#                         baset_heating={
#                             "paved": 0.1,
#                             "bldgs": 0.1,
#                             "dectr": 0.1,
#                             "evetr": 0.1,
#                             "grass": 0.1,
#                             "bsoil": 0.1,
#                             "water": 0.1,
#                         },
#                         ah_min=0.1,
#                         ah_slope_cooling=0.1,
#                         ah_slope_heating=0.1,
#                         ahprof_24hr=HourlyProfile(
#                             working_day={str(i): 0.1 for i in range(1, 25)},
#                             holiday={str(i): 0.1 for i in range(1, 25)},
#                         ),
#                     ),
#                     co2=CO2Params(
#                         co2pointsource=0.1,
#                         ef_umolco2perj=0.1,
#                         enef_v_jkm=0.1,
#                         fcef_v_kgkm=0.1,
#                         frfossilfuel_heat=0.1,
#                         frfossilfuel_nonheat=0.1,
#                         maxfcmetab=0.1,
#                         maxqfmetab=0.1,
#                         min_res_bioco2=0.1,
#                         minfcmetab=0.1,
#                         minqfmetab=0.1,
#                         trafficrate=0.1,
#                         trafficunits=0.1,
#                         traffprof_24hr=HourlyProfile(
#                             working_day={str(i): 0.1 for i in range(1, 25)},
#                             holiday={str(i): 0.1 for i in range(1, 25)},
#                         ),
#                     ),
#                 ),
#                 snow=SnowParams(
#                     crwmax=0.1,
#                     crwmin=0.1,
#                     narp_emis_snow=0.1,
#                     preciplimit=0.1,
#                     preciplimitalb=0.1,
#                     snowalbmax=0.1,
#                     snowalbmin=0.1,
#                     snowdensmin=0.1,
#                     snowlimbldg=0.1,
#                     snowlimpaved=0.1,
#                     snowpacklimit=0.1,
#                     snowprof_24hr=HourlyProfile(
#                         working_day={str(i): 0.1 for i in range(1, 25)},
#                         holiday={str(i): 0.1 for i in range(1, 25)},
#                     ),
#                     tau_a=0.1,
#                     tau_r=0.1,
#                     tempmeltfact=0.1,
#                     radmeltfact=0.1,
#                 ),
#                 land_cover=LandCover(
#                     paved=SurfaceProperties(
#                         sfr=0.1,
#                         emis=0.1,
#                         alb=0.1,
#                         chanohm=0.1,
#                         cpanohm=0.1,
#                         kkanohm=0.1,
#                         ohm_threshsw=0.1,
#                         ohm_threshwd=0.1,
#                         soildepth=0.1,
#                         soilstore=0.1,
#                         soilstorecap=0.1,
#                         statelimit=0.1,
#                         sathydraulicconduct=0.1,
#                         storedrainprm=StorageDrainParams(
#                             store_min=0.0,
#                             store_max=10.0,
#                             store_cap=5.0,
#                             drain_eq=1,
#                             drain_coef_1=0.5,
#                             drain_coef_2=0.5,
#                         ),
#                     ),
#                     bldgs=SurfaceProperties(
#                         sfr=0.1,
#                         emis=0.1,
#                         alb=0.1,
#                         chanohm=0.1,
#                         cpanohm=0.1,
#                         kkanohm=0.1,
#                         ohm_threshsw=0.1,
#                         ohm_threshwd=0.1,
#                         soildepth=0.1,
#                         soilstore=0.1,
#                         soilstorecap=0.1,
#                         statelimit=0.1,
#                         sathydraulicconduct=0.1,
#                         storedrainprm=StorageDrainParams(
#                             store_min=0.0,
#                             store_max=10.0,
#                             store_cap=5.0,
#                             drain_eq=1,
#                             drain_coef_1=0.5,
#                             drain_coef_2=0.5,
#                         ),
#                     ),
#                     dectr=DectrProperties(
#                         sfr=0.1,
#                         emis=0.1,
#                         alb=0.1,
#                         soildepth=0.1,
#                         soilstore=0.1,
#                         soilstorecap=0.1,
#                         statelimit=0.1,
#                         sathydraulicconduct=0.1,
#                         storedrainprm=StorageDrainParams(
#                             store_min=0.0,
#                             store_max=10.0,
#                             store_cap=5.0,
#                             drain_eq=1,
#                             drain_coef_1=0.5,
#                             drain_coef_2=0.5,
#                         ),
#                         faidectree=0.1,
#                         dectreeh=0.1,
#                         pormin_dec=0.1,
#                         pormax_dec=0.1,
#                         capmax_dec=0.1,
#                         capmin_dec=0.1,
#                         alb_min=0.1,
#                         alb_max=0.9,
#                         beta_bioco2=0.1,
#                         beta_enh_bioco2=0.1,
#                         alpha_bioco2=0.1,
#                         alpha_enh_bioco2=0.1,
#                         resp_a=0.1,
#                         resp_b=0.1,
#                         theta_bioco2=0.1,
#                         conductance=Conductance(
#                             g_max=0.1,
#                             g_k=0.1,
#                             g_q_base=0.1,
#                             g_q_shape=0.1,
#                             g_t=0.1,
#                             g_sm=0.1,
#                             kmax=0.1,
#                             gsmodel=1,
#                             maxconductance=0.1,
#                             s1=0.1,
#                             s2=0.1,
#                         ),
#                         lai=LAIParams(
#                             baset=0.1,
#                             gddfull=0.1,
#                             basete=0.1,
#                             sddfull=0.1,
#                             laimin=0.1,
#                             laimax=0.1,
#                             laipower=0.1,
#                             laitype=1,
#                         ),
#                         ie_a=0.1,
#                         ie_m=0.1,
#                     ),
#                     evetr=EvetrProperties(
#                         sfr=0.1,
#                         emis=0.1,
#                         alb=0.1,
#                         soildepth=0.1,
#                         soilstore=0.1,
#                         soilstorecap=0.1,
#                         statelimit=0.1,
#                         sathydraulicconduct=0.1,
#                         storedrainprm=StorageDrainParams(
#                             store_min=0.0,
#                             store_max=10.0,
#                             store_cap=5.0,
#                             drain_eq=1,
#                             drain_coef_1=0.5,
#                             drain_coef_2=0.5,
#                         ),
#                         faievetree=0.1,
#                         evetreeh=0.1,
#                         alb_min=0.1,
#                         alb_max=0.9,
#                         beta_bioco2=0.1,
#                         beta_enh_bioco2=0.1,
#                         alpha_bioco2=0.1,
#                         alpha_enh_bioco2=0.1,
#                         resp_a=0.1,
#                         resp_b=0.1,
#                         theta_bioco2=0.1,
#                         conductance=Conductance(
#                             g_max=0.1,
#                             g_k=0.1,
#                             g_q_base=0.1,
#                             g_q_shape=0.1,
#                             g_t=0.1,
#                             g_sm=0.1,
#                             kmax=0.1,
#                             gsmodel=1,
#                             maxconductance=0.1,
#                             s1=0.1,
#                             s2=0.1,
#                         ),
#                         lai=LAIParams(
#                             baset=0.1,
#                             gddfull=0.1,
#                             basete=0.1,
#                             sddfull=0.1,
#                             laimin=0.1,
#                             laimax=0.1,
#                             laipower=0.1,
#                             laitype=1,
#                         ),
#                         ie_a=0.1,
#                         ie_m=0.1,
#                     ),
#                     grass=VegetatedSurfaceProperties(
#                         sfr=0.1,
#                         emis=0.1,
#                         alb=0.1,
#                         soildepth=0.1,
#                         soilstore=0.1,
#                         soilstorecap=0.1,
#                         statelimit=0.1,
#                         sathydraulicconduct=0.1,
#                         storedrainprm=StorageDrainParams(
#                             store_min=0.0,
#                             store_max=10.0,
#                             store_cap=5.0,
#                             drain_eq=1,
#                             drain_coef_1=0.5,
#                             drain_coef_2=0.5,
#                         ),
#                         alb_min=0.1,
#                         alb_max=0.9,
#                         beta_bioco2=0.1,
#                         beta_enh_bioco2=0.1,
#                         alpha_bioco2=0.1,
#                         alpha_enh_bioco2=0.1,
#                         resp_a=0.1,
#                         resp_b=0.1,
#                         theta_bioco2=0.1,
#                         conductance=Conductance(
#                             g_max=0.1,
#                             g_k=0.1,
#                             g_q_base=0.1,
#                             g_q_shape=0.1,
#                             g_t=0.1,
#                             g_sm=0.1,
#                             kmax=0.1,
#                             gsmodel=1,
#                             maxconductance=0.1,
#                             s1=0.1,
#                             s2=0.1,
#                         ),
#                         lai=LAIParams(
#                             baset=0.1,
#                             gddfull=0.1,
#                             basete=0.1,
#                             sddfull=0.1,
#                             laimin=0.1,
#                             laimax=0.1,
#                             laipower=0.1,
#                             laitype=1,
#                         ),
#                         ie_a=0.1,
#                         ie_m=0.1,
#                     ),
#                     bsoil=SurfaceProperties(
#                         sfr=0.1,
#                         emis=0.1,
#                         alb=0.1,
#                         chanohm=0.1,
#                         cpanohm=0.1,
#                         kkanohm=0.1,
#                         ohm_threshsw=0.1,
#                         ohm_threshwd=0.1,
#                         soildepth=0.1,
#                         soilstore=0.1,
#                         soilstorecap=0.1,
#                         statelimit=0.1,
#                         sathydraulicconduct=0.1,
#                         storedrainprm=StorageDrainParams(
#                             store_min=0.0,
#                             store_max=10.0,
#                             store_cap=5.0,
#                             drain_eq=1,
#                             drain_coef_1=0.5,
#                             drain_coef_2=0.5,
#                         ),
#                     ),
#                     water=WaterProperties(
#                         sfr=0.1,
#                         emis=0.1,
#                         alb=0.1,
#                         soildepth=0.1,
#                         soilstore=0.1,
#                         soilstorecap=0.1,
#                         statelimit=0.1,
#                         sathydraulicconduct=0.1,
#                         storedrainprm=StorageDrainParams(
#                             store_min=0.0,
#                             store_max=10.0,
#                             store_cap=5.0,
#                             drain_eq=1,
#                             drain_coef_1=0.5,
#                             drain_coef_2=0.5,
#                         ),
#                         flowchange=0.1,
#                         irrfracwater=0.1,
#                     ),
#                 ),
#                 initial_states={
#                     "paved": InitialState(
#                         state=0.1,
#                         soilstore=0.1,
#                         snowfrac=0.1,
#                         snowpack=0.1,
#                         snowwater=0.1,
#                         snowalb=0.1,
#                         snowdens=0.1,
#                     ),
#                     "bldgs": InitialState(
#                         state=0.1,
#                         soilstore=0.1,
#                         snowfrac=0.1,
#                         snowpack=0.1,
#                         snowwater=0.1,
#                         snowalb=0.1,
#                         snowdens=0.1,
#                     ),
#                     "dectr": InitialState(
#                         state=0.1,
#                         soilstore=0.1,
#                         snowfrac=0.1,
#                         snowpack=0.1,
#                         snowwater=0.1,
#                         snowalb=0.1,
#                         snowdens=0.1,
#                         decidcap_id=0.1,
#                         porosity_id=0.1,
#                         alb_id=0.1,
#                     ),
#                     "evetr": InitialState(
#                         state=0.1,
#                         soilstore=0.1,
#                         snowfrac=0.1,
#                         snowpack=0.1,
#                         snowwater=0.1,
#                         snowalb=0.1,
#                         snowdens=0.1,
#                         decidcap_id=0.1,
#                         porosity_id=0.1,
#                         alb_id=0.1,
#                     ),
#                     "grass": InitialState(
#                         state=0.1,
#                         soilstore=0.1,
#                         snowfrac=0.1,
#                         snowpack=0.1,
#                         snowwater=0.1,
#                         snowalb=0.1,
#                         snowdens=0.1,
#                         decidcap_id=0.1,
#                         porosity_id=0.1,
#                         alb_id=0.1,
#                     ),
#                     "bsoil": InitialState(
#                         state=0.1,
#                         soilstore=0.1,
#                         snowfrac=0.1,
#                         snowpack=0.1,
#                         snowwater=0.1,
#                         snowalb=0.1,
#                         snowdens=0.1,
#                     ),
#                     "water": InitialState(
#                         state=0.1,
#                         soilstore=0.1,
#                         snowfrac=0.1,
#                         snowpack=0.1,
#                         snowwater=0.1,
#                         snowalb=0.1,
#                         snowdens=0.1,
#                     ),
#                 },
#             ),
#         )
#     ],
# )


if __name__ == "__main__":
    # test the sample config
    # Load YAML config
    with open('./config-suews.yml', 'r') as file:
        yaml_config = yaml.safe_load(file)

    # Create SUEWSConfig object
    suews_config = SUEWSConfig(**yaml_config[0])

    print(r'testing suews_config done!')
    # test converting the sample `config-suews.yml` to dataframe



