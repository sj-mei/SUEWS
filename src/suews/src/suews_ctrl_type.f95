MODULE SUEWS_DEF_DTS
    USE allocateArray, ONLY: nsurf, nvegsurf, ndepth
    ! ********** SUEWS_parameters schema (basic) **********
    TYPE, PUBLIC :: METHOD_PRM
       INTEGER :: DiagMethod ! Defines how near surface diagnostics are calculated
       INTEGER :: EmissionsMethod ! method to calculate anthropogenic heat [-]
       INTEGER :: RoughLenHeatMethod ! method to calculate heat roughness length [-]
       INTEGER :: RoughLenMomMethod ! Determines how aerodynamic roughness length (z0m) and zero displacement height (zdm) are calculated [-]
       INTEGER :: FAIMethod !Determines how FAI is calculated [-]
       INTEGER :: SMDMethod ! Determines method for calculating soil moisture deficit [-]
       INTEGER :: WaterUseMethod ! Defines how external water use is calculated[-]
       INTEGER :: NetRadiationMethod ! method for calculation of radiation fluxes [-]
       INTEGER :: StabilityMethod ! method to calculate atmospheric stability [-]
       INTEGER :: StorageHeatMethod ! !Determines method for calculating storage heat flux ΔQS [-]
       INTEGER :: Diagnose ! flag for printing diagnostic info during runtime [N/A]C
       INTEGER :: SnowUse !
       LOGICAL :: use_sw_direct_albedo !boolean, Specify ground and roof albedos separately for direct solar radiation [-]
       INTEGER :: ohmIncQF ! Determines whether the storage heat flux calculation uses Q* or ( Q* +QF) [-]
    END TYPE METHOD_PRM
 
    TYPE, PUBLIC :: SURF_STORE_PRM
       REAL(KIND(1D0)) :: store_min
       REAL(KIND(1D0)) :: store_max
       REAL(KIND(1D0)) :: store_cap ! this should be abundant ? (current storage capacity) [mm]
       INTEGER :: drain_eq
       REAL(KIND(1D0)) :: drain_coef_1
       REAL(KIND(1D0)) :: drain_coef_2
    END TYPE SURF_STORE_PRM
 
    TYPE, PUBLIC :: WATER_DIST_PRM
       REAL(KIND(1D0)) :: to_paved
       REAL(KIND(1D0)) :: to_bldg
       REAL(KIND(1D0)) :: to_evetr
       REAL(KIND(1D0)) :: to_dectr
       REAL(KIND(1D0)) :: to_grass
       REAL(KIND(1D0)) :: to_bsoil
       REAL(KIND(1D0)) :: to_water
       REAL(KIND(1D0)) :: to_soilstore
    END TYPE WATER_DIST_PRM
 
    TYPE, PUBLIC :: bioCO2_PRM
       REAL(KIND(1D0)) :: beta_bioco2 ! The light-saturated gross photosynthesis of the canopy [umol m-2 s-1 ]
       REAL(KIND(1D0)) :: beta_enh_bioco2 ! Part of the beta coefficient related to the fraction of vegetation [umol m-2 s-1 ]
       REAL(KIND(1D0)) :: alpha_bioco2 ! The mean apparent ecosystem quantum. Represents the initial slope of the light-response curve [-]
       REAL(KIND(1D0)) :: alpha_enh_bioco2 ! Part of the alpha coefficient related to the fraction of vegetation [-]
       REAL(KIND(1D0)) :: resp_a !Respiration coefficient a
       REAL(KIND(1D0)) :: resp_b !Respiration coefficient b - related to air temperature dependency
       REAL(KIND(1D0)) :: theta_bioco2 ! The convexity of the curve at light saturation [-]
       REAL(KIND(1D0)) :: min_res_bioCO2 ! Minimum soil respiration rate (for cold-temperature limit) [umol m-2 s-1]
    END TYPE bioCO2_PRM
 
    TYPE, PUBLIC :: CONDUCTANCE_PRM
       REAL(KIND(1D0)) :: g_max !Fitted parameters related to surface res. calculations
       REAL(KIND(1D0)) :: g_k
       REAL(KIND(1D0)) :: g_q_base
       REAL(KIND(1D0)) :: g_q_shape
       REAL(KIND(1D0)) :: g_t
       REAL(KIND(1D0)) :: g_sm ! Fitted parameters related to surface res. calculations
       REAL(KIND(1D0)) :: kmax ! annual maximum hourly solar radiation [W m-2]
       INTEGER :: gsmodel ! choice of gs parameterisation (1 = Ja11, 2 = Wa16) [-]
       REAL(KIND(1D0)) :: s1 ! a parameter related to soil moisture dependence [-]
       REAL(KIND(1D0)) :: s2 ! a parameter related to soil moisture dependence [mm]
       REAL(KIND(1D0)) :: TH ! upper air temperature limit [degC]
       REAL(KIND(1D0)) :: TL ! lower air temperature limit [degC]
    END TYPE CONDUCTANCE_PRM
 
    TYPE, PUBLIC :: LAI_PRM ! do we need `lai_id` here?
       REAL(KIND(1D0)) :: baset !Base Temperature for initiating growing degree days (GDD) for leaf growth [degC]
       REAL(KIND(1D0)) :: gddfull !the growing degree days (GDD) needed for full capacity of the leaf area index [degC]
       REAL(KIND(1D0)) :: basete !Base temperature for initiating sensesance degree days (SDD) for leaf off [degC]
       REAL(KIND(1D0)) :: sddfull !the sensesence degree days (SDD) needed to initiate leaf off [degC]
       REAL(KIND(1D0)) :: laimin !leaf-off wintertime value [m2 m-2]
       REAL(KIND(1D0)) :: laimax !full leaf-on summertime value [m2 m-2]
       REAL(KIND(1D0)), DIMENSION(4) :: laipower ! parameters required by LAI calculation.
       INTEGER :: laitype ! LAI calculation choice.
    END TYPE LAI_PRM
 
    TYPE, PUBLIC :: OHM_COEF_LC
       REAL(KIND(1D0)) :: summer_dry
       REAL(KIND(1D0)) :: summer_wet
       REAL(KIND(1D0)) :: winter_dry
       REAL(KIND(1D0)) :: winter_wet
    END TYPE OHM_COEF_LC
 
    TYPE, PUBLIC :: OHM_PRM
       REAL(KIND(1D0)) :: chanohm ! Bulk transfer coefficient for this surface to use in AnOHM [J m-3 K-1]
       REAL(KIND(1D0)) :: cpanohm ! Volumetric heat capacity for this surface to use in AnOHM  [J m-3 K-1]
       REAL(KIND(1D0)) :: kkanohm ! Thermal conductivity for this surface to use in AnOHM [W m-1 K-1]
       REAL(KIND(1D0)) :: ohm_threshsw ! Temperature threshold determining whether summer/winter OHM coefficients are applied [degC]
       REAL(KIND(1D0)) :: ohm_threshwd ! Soil moisture threshold determining whether wet/dry OHM coefficients are applied [degC]
       TYPE(OHM_COEF_LC), DIMENSION(3) :: ohm_coef_lc
    END TYPE OHM_PRM
 
    TYPE, PUBLIC :: SOIL_PRM
       REAL(KIND(1D0)) :: soildepth
       REAL(KIND(1D0)) :: soilstorecap
       REAL(KIND(1D0)) :: sathydraulicconduct ! !Hydraulic conductivity for saturated soil [mm s-1]
    END TYPE SOIL_PRM
 
    TYPE, PUBLIC :: anthroHEAT_PRM
       REAL(KIND(1D0)) :: qf0_beu_working ! Fraction of base value coming from buildings (working day) [-]
       REAL(KIND(1D0)) :: qf0_beu_holiday ! Fraction of base value coming from buildings (holiday) [-]
       REAL(KIND(1D0)) :: qf_a_working ! Base value for QF (working day) [W m-2]
       REAL(KIND(1D0)) :: qf_a_holiday ! Base value for QF (holiday) [W m-2]
       REAL(KIND(1D0)) :: qf_b_working ! Parameter related to heating degree days (working day) [W m-2 K-1 (Cap ha-1 )-1]
       REAL(KIND(1D0)) :: qf_b_holiday ! Parameter related to heating degree days (holiday) [W m-2 K-1 (Cap ha-1 )-1]
       REAL(KIND(1D0)) :: qf_c_working ! Parameter related to cooling degree days (working day) [W m-2 K-1 (Cap ha-1 )-1]
       REAL(KIND(1D0)) :: qf_c_holiday ! Parameter related to cooling degree days (holiday) [W m-2 K-1 (Cap ha-1 )-1]
       REAL(KIND(1D0)) :: baset_cooling_working ! Base temperature for cooling degree days (working day) [degC]
       REAL(KIND(1D0)) :: baset_cooling_holiday ! Base temperature for cooling degree days (holiday) [degC]
       REAL(KIND(1D0)) :: baset_heating_working ! Base temperature for heating degree days (working day) [degC]
       REAL(KIND(1D0)) :: baset_heating_holiday ! Base temperature for heating degree days (holiday) [degC]
       REAL(KIND(1D0)) :: ah_min_working ! minimum QF values (working day) [W m-2]
       REAL(KIND(1D0)) :: ah_min_holiday ! minimum QF values (holiday) [W m-2]
       REAL(KIND(1D0)) :: ah_slope_cooling_working ! cooling slope for the anthropogenic heat flux calculation (working day) [W m-2 K-1]
       REAL(KIND(1D0)) :: ah_slope_cooling_holiday ! cooling slope for the anthropogenic heat flux calculation (holiday) [W m-2 K-1]
       REAL(KIND(1D0)) :: ah_slope_heating_working ! heating slope for the anthropogenic heat flux calculation (working day) [W m-2 K-1]
       REAL(KIND(1D0)) :: ah_slope_heating_holiday ! heating slope for the anthropogenic heat flux calculation (holiday) [W m-2 K-1]
       REAL(KIND(1D0)), DIMENSION(0:23) :: ahprof_24hr_working ! Hourly profile values used in energy use calculation (working day) [-]
       REAL(KIND(1D0)), DIMENSION(0:23) :: ahprof_24hr_holiday ! Hourly profile values used in energy use calculation (holiday) [-]
       REAL(KIND(1D0)) :: popdensdaytime_working ! Daytime population density [people ha-1] (working day)
       REAL(KIND(1D0)) :: popdensdaytime_holiday ! Daytime population density [people ha-1] (holiday)
       REAL(KIND(1D0)) :: popdensnighttime
       REAL(KIND(1D0)), DIMENSION(0:23) :: popprof_24hr_working !Hourly profile values used in dynamic population estimation[-] (working day)
       REAL(KIND(1D0)), DIMENSION(0:23) :: popprof_24hr_holiday !Hourly profile values used in dynamic population estimation[-] (holiday)
    END TYPE anthroHEAT_PRM
 
    TYPE, PUBLIC :: IRRIG_daywater
       REAL(KIND(1D0)) :: monday_flag ! Irrigation flag: 1 for on and 0 for off.
       REAL(KIND(1D0)) :: monday_percent ! Fraction of properties using irrigation for each day of a week.
       REAL(KIND(1D0)) :: tuesday_flag
       REAL(KIND(1D0)) :: tuesday_percent
       REAL(KIND(1D0)) :: wednesday_flag
       REAL(KIND(1D0)) :: wednesday_percent
       REAL(KIND(1D0)) :: thursday_flag
       REAL(KIND(1D0)) :: thursday_percent
       REAL(KIND(1D0)) :: friday_flag
       REAL(KIND(1D0)) :: friday_percent
       REAL(KIND(1D0)) :: saturday_flag
       REAL(KIND(1D0)) :: saturday_percent
       REAL(KIND(1D0)) :: sunday_flag
       REAL(KIND(1D0)) :: sunday_percent
    END TYPE IRRIG_daywater
 
    TYPE, PUBLIC :: IRRIGATION_PRM ! used in irrigation
       REAL(KIND(1D0)) :: h_maintain ! ponding water depth to maintain [mm]
       REAL(KIND(1D0)) :: faut ! Fraction of irrigated area using automatic irrigation [-]
       REAL(KIND(1D0)), DIMENSION(3) :: ie_a ! Coefficient for automatic irrigation model
       REAL(KIND(1D0)), DIMENSION(3) :: ie_m ! Coefficient for automatic irrigation model
       INTEGER :: ie_start ! ending time of water use [DOY]
       INTEGER :: ie_end ! starting time of water use [DOY]
       REAL(KIND(1D0)) :: internalwateruse_h ! Internal water use [mm h-1]
       TYPE(IRRIG_daywater) :: irr_daywater
       REAL(KIND(1D0)), DIMENSION(0:23) :: wuprofa_24hr_working ! Hourly profile values used in automatic irrigation[-]
       REAL(KIND(1D0)), DIMENSION(0:23) :: wuprofa_24hr_holiday
       REAL(KIND(1D0)), DIMENSION(0:23) :: wuprofm_24hr_working ! Hourly profile values used in manual irrigation[-]
       REAL(KIND(1D0)), DIMENSION(0:23) :: wuprofm_24hr_holiday
    END TYPE IRRIGATION_PRM
 
    TYPE, PUBLIC :: anthroEMIS_PRM
       INTEGER :: startdls ! start of daylight saving  [DOY]
       INTEGER :: enddls ! end of daylight saving [DOY]
       TYPE(anthroHEAT_PRM) :: anthroheat
       REAL(KIND(1D0)) :: EF_umolCO2perJ
       REAL(KIND(1D0)) :: EnEF_v_Jkm
       REAL(KIND(1D0)) :: FrFossilFuel_Heat
       REAL(KIND(1D0)) :: FrFossilFuel_NonHeat
       REAL(KIND(1D0)), DIMENSION(2) :: FcEF_v_kgkm
       REAL(KIND(1D0)), DIMENSION(0:23) :: HumActivity_24hr_working
       REAL(KIND(1D0)), DIMENSION(0:23) :: HumActivity_24hr_holiday
       REAL(KIND(1D0)) :: MaxFCMetab
       REAL(KIND(1D0)) :: MaxQFMetab
       REAL(KIND(1D0)) :: MinFCMetab
       REAL(KIND(1D0)) :: MinQFMetab
       REAL(KIND(1D0)) :: TrafficRate_working
       REAL(KIND(1D0)) :: TrafficRate_holiday
       REAL(KIND(1D0)) :: TrafficUnits
       REAL(KIND(1D0)), DIMENSION(0:23) :: TraffProf_24hr_working
       REAL(KIND(1D0)), DIMENSION(0:23) :: TraffProf_24hr_holiday
    END TYPE anthroEMIS_PRM
 
    TYPE, PUBLIC :: SNOW_PRM
       REAL(KIND(1D0)) :: crwmax ! maximum water holding capacity of snow [mm]
       REAL(KIND(1D0)) :: crwmin ! minimum water holding capacity of snow [mm]
       REAL(KIND(1D0)) :: narp_emis_snow ! snow emissivity in NARP model [-]
       REAL(KIND(1D0)) :: preciplimit ! temperature limit when precipitation falls as snow [degC]
       REAL(KIND(1D0)) :: preciplimitalb ! Limit for hourly precipitation when the ground is fully covered with snow [mm]
       REAL(KIND(1D0)) :: snowalbmax ! effective surface albedo (middle of the day value) for summertime [-]
       REAL(KIND(1D0)) :: snowalbmin ! effective surface albedo (middle of the day value) for wintertime (not including snow) [-]
       REAL(KIND(1D0)) :: snowdensmax ! maximum snow density [kg m-3]
       REAL(KIND(1D0)) :: snowdensmin ! fresh snow density [kg m-3]
       REAL(KIND(1D0)) :: snowlimbldg ! Limit of the snow water equivalent for snow removal from building roofs [mm]
       REAL(KIND(1D0)) :: snowlimpaved ! limit of the snow water equivalent for snow removal from roads[mm]
       REAL(KIND(1D0)), DIMENSION(nsurf) :: snowpacklimit ! Limit for the snow water equivalent when snow cover starts to be patchy [mm]
       REAL(KIND(1D0)), DIMENSION(0:23) :: snowprof_24hr_working ! Hourly profile values used in snow clearing of working day [-]
       REAL(KIND(1D0)), DIMENSION(0:23) :: snowprof_24hr_holiday ! Hourly profile values used in snow clearing of holiday [-]
       REAL(KIND(1D0)) :: tau_a ! time constant for snow albedo aging in cold snow [-]
       REAL(KIND(1D0)) :: tau_f ! time constant for snow albedo aging in melting snow [-]
       REAL(KIND(1D0)) :: tau_r ! time constant for snow density ageing [-]
       REAL(KIND(1D0)) :: tempmeltfact ! hourly temperature melt factor of snow [mm K-1 h-1]
       REAL(KIND(1D0)) :: radmeltfact ! hourly radiation melt factor of snow [mm W-1 h-1]
    END TYPE SNOW_PRM
 
    TYPE, PUBLIC :: SPARTACUS_PRM
       REAL(KIND(1D0)) :: air_ext_lw
       REAL(KIND(1D0)) :: air_ext_sw
       REAL(KIND(1D0)) :: air_ssa_lw
       REAL(KIND(1D0)) :: air_ssa_sw
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: height
       REAL(KIND(1D0)) :: ground_albedo_dir_mult_fact
       INTEGER :: n_stream_lw_urban ! LW streams per hemisphere [-]
       INTEGER :: n_stream_sw_urban ! shortwave diffuse streams per hemisphere [-]
       INTEGER :: n_vegetation_region_urban ! Number of regions used to describe vegetation [-]
       REAL(KIND(1D0)) :: sw_dn_direct_frac
       REAL(KIND(1D0)) :: use_sw_direct_albedo
       REAL(KIND(1D0)) :: veg_contact_fraction_const
       REAL(KIND(1D0)) :: veg_fsd_const
       REAL(KIND(1D0)) :: veg_ssa_lw
       REAL(KIND(1D0)) :: veg_ssa_sw
    END TYPE SPARTACUS_PRM
 
    TYPE, PUBLIC :: SPARTACUS_LAYER_PRM
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: building_frac ! building fraction [-]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: building_scale ! diameter of buildings [[m]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: veg_frac ! vegetation fraction [-]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: veg_scale ! scale of tree crowns [m]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: alb_roof ! albedo of roof [-]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: emis_roof ! emissivity of roof [-]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: alb_wall ! albedo of wall [-]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: emis_wall ! emissivity of wall [-]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: roof_albedo_dir_mult_fact !Ratio of the direct and diffuse albedo of the roof[-]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: wall_specular_frac ! Fraction of wall reflection that is specular [-]
    END TYPE SPARTACUS_LAYER_PRM
 
    ! ********** SUEWS_parameters schema (derived) **********
    TYPE, PUBLIC :: SITE_PRM
       REAL(KIND(1D0)) :: lat !latitude [deg]
       REAL(KIND(1D0)) :: lon !longitude [deg]
       REAL(KIND(1D0)) :: alt ! solar altitude [deg]
       INTEGER :: gridiv ! grid id [-]
       REAL(KIND(1D0)) :: timezone ! time zone, for site relative to UTC (east is positive) [h]
       REAL(KIND(1D0)) :: surfacearea ! area of the grid [ha]
       REAL(KIND(1D0)) :: z ! measurement height [m]
       REAL(KIND(1D0)) :: z0m_in ! roughness length for momentum [m]
       REAL(KIND(1D0)) :: zdm_in ! zero-plane displacement [m]
       REAL(KIND(1D0)) :: pipecapacity ! capacity of pipes to transfer water [mm]
       REAL(KIND(1D0)) :: runofftowater ! fraction of above-ground runoff flowing to water surface during flooding [-]
       REAL(KIND(1D0)) :: narp_trans_site ! atmospheric transmissivity for NARP [-]
       REAL(KIND(1D0)) :: CO2PointSource ! CO2 emission factor [kg km-1]
       REAL(KIND(1D0)) :: flowchange ! Difference in input and output flows for water surface
    END TYPE SITE_PRM
 
    TYPE, PUBLIC :: LUMPS_PRM
       REAL(KIND(1D0)) :: raincover ! limit when surface totally covered with water for LUMPS [mm]
       REAL(KIND(1D0)) :: rainmaxres ! maximum water bucket reservoir. Used for LUMPS surface wetness control. [mm]
       REAL(KIND(1D0)) :: drainrt ! Drainage rate of the water bucket [mm hr-1]
       INTEGER :: veg_type ! Defines how vegetation is calculated for LUMPS [-]
    END TYPE LUMPS_PRM
 
    TYPE, PUBLIC :: EHC_PRM
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: soil_storecap_roof ! Capacity of soil store for roof [mm]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: soil_storecap_wall ! Capacity of soil store for wall [mm]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: state_limit_roof ! Limit for state_id of roof [mm]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: state_limit_wall ! Limit for state_id of wall [mm]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: wet_thresh_roof ! wetness threshold  of roof [mm]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: wet_thresh_wall ! wetness threshold  of wall [mm]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: tin_roof ! indoor temperature for roof [degC]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: tin_wall ! indoor temperature for wall [degC]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: tin_surf ! deep bottom temperature for each surface [degC]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: k_roof ! thermal conductivity of roof [W m-1 K]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: k_wall ! thermal conductivity of wall [W m-1 K]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: k_surf ! thermal conductivity of v [W m-1 K]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: cp_roof ! Heat capacity of roof [J m-3 K-1]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: cp_wall ! Heat capacity of wall [J m-3 K-1]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: cp_surf ! Heat capacity of each surface [J m-3 K-1]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: dz_roof ! thickness of each layer in roof [m]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: dz_wall ! thickness of each layer in wall [m]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: dz_surf ! thickness of each layer in surface [m]
    END TYPE EHC_PRM
 
    TYPE, PUBLIC :: LC_PAVED_PRM
       REAL(KIND(1D0)) :: sfr
       REAL(KIND(1D0)) :: emis
       TYPE(OHM_PRM) :: ohm
       TYPE(SOIL_PRM) :: soil
       REAL(KIND(1D0)) :: state
       REAL(KIND(1D0)) :: statelimit
       REAL(KIND(1D0)) :: irrfracpaved
       REAL(KIND(1D0)) :: wetthresh
       TYPE(WATER_DIST_PRM) :: waterdist
    END TYPE LC_PAVED_PRM
 
    TYPE, PUBLIC :: LC_BLDG_PRM
       REAL(KIND(1D0)) :: sfr
       REAL(KIND(1D0)) :: faibldg
       REAL(KIND(1D0)) :: bldgh
       REAL(KIND(1D0)) :: emis
       TYPE(OHM_PRM) :: ohm
       TYPE(SOIL_PRM) :: soil
       REAL(KIND(1D0)) :: state
       REAL(KIND(1D0)) :: statelimit
       REAL(KIND(1D0)) :: irrfracbldgs
       REAL(KIND(1D0)) :: wetthresh
       TYPE(WATER_DIST_PRM) :: waterdist
    END TYPE LC_BLDG_PRM
 
    TYPE, PUBLIC :: LC_DECTR_PRM
       REAL(KIND(1D0)) :: sfr
       REAL(KIND(1D0)) :: emis
       REAL(KIND(1D0)) :: faidectree
       REAL(KIND(1D0)) :: dectreeh
       REAL(KIND(1D0)) :: pormin_dec ! absent for evergreen trees ??
       REAL(KIND(1D0)) :: pormax_dec
       REAL(KIND(1D0)) :: alb_min
       REAL(KIND(1D0)) :: alb_max
       TYPE(OHM_PRM) :: ohm
       TYPE(SOIL_PRM) :: soil
       REAL(KIND(1D0)) :: statelimit ! ******* dummy variable *******
       REAL(KIND(1D0)) :: capmax_dec ! Maximum water storage capacity for upper surfaces (i.e. canopy) (absent for evergreen trees ??)
       REAL(KIND(1D0)) :: capmin_dec ! Minimum water storage capacity for upper surfaces (i.e. canopy).
       REAL(KIND(1D0)) :: irrfracdectr
       REAL(KIND(1D0)) :: wetthresh
       TYPE(bioCO2_PRM) :: bioco2
       REAL(KIND(1D0)) :: maxconductance ! the maximum conductance of each vegetation or surface type. [mm s-1]
       ! TYPE(CONDUCTANCE_PRM) :: conductance
       TYPE(LAI_PRM) :: lai
       TYPE(WATER_DIST_PRM) :: waterdist
    END TYPE LC_DECTR_PRM
 
    TYPE, PUBLIC :: LC_EVETR_PRM
       REAL(KIND(1D0)) :: sfr !surface cover fraction[-]
       REAL(KIND(1D0)) :: emis !Effective surface emissivity[-]
       REAL(KIND(1D0)) :: faievetree !frontal area index for evergreen tree [-]
       REAL(KIND(1D0)) :: evetreeh !height of evergreen tree [m]
       REAL(KIND(1D0)) :: alb_min !minimum albedo for evergreen tree [-]
       REAL(KIND(1D0)) :: alb_max !maximum albedo for evergreen tree [-]
       TYPE(OHM_PRM) :: ohm
       TYPE(SOIL_PRM) :: soil
       REAL(KIND(1D0)) :: statelimit ! upper limit to the surface state [mm]
       REAL(KIND(1D0)) :: irrfracevetr
       REAL(KIND(1D0)) :: wetthresh !surface wetness threshold [mm], When State > WetThresh, RS=0 limit in SUEWS_evap [mm]
       TYPE(bioCO2_PRM) :: bioco2
       ! TYPE(CONDUCTANCE_PRM) :: conductance
       REAL(KIND(1D0)) :: maxconductance ! the maximum conductance of each vegetation or surface type. [mm s-1]
       TYPE(LAI_PRM) :: lai
       TYPE(WATER_DIST_PRM) :: waterdist !Fraction of water redistribution [-]
    END TYPE LC_EVETR_PRM
 
    TYPE, PUBLIC :: LC_GRASS_PRM
       REAL(KIND(1D0)) :: sfr
       REAL(KIND(1D0)) :: emis
       REAL(KIND(1D0)) :: alb_min
       REAL(KIND(1D0)) :: alb_max
       TYPE(OHM_PRM) :: ohm
       TYPE(SOIL_PRM) :: soil
       REAL(KIND(1D0)) :: statelimit
       REAL(KIND(1D0)) :: irrfracgrass
       REAL(KIND(1D0)) :: wetthresh
       TYPE(bioCO2_PRM) :: bioco2
       ! TYPE(CONDUCTANCE_PRM) :: conductance
       REAL(KIND(1D0)) :: maxconductance ! the maximum conductance of each vegetation or surface type. [mm s-1]
       TYPE(LAI_PRM) :: lai
       TYPE(WATER_DIST_PRM) :: waterdist
    END TYPE LC_GRASS_PRM
 
    TYPE, PUBLIC :: LC_BSOIL_PRM
       REAL(KIND(1D0)) :: sfr
       REAL(KIND(1D0)) :: emis
       TYPE(OHM_PRM) :: ohm
       TYPE(SOIL_PRM) :: soil
       REAL(KIND(1D0)) :: statelimit
       REAL(KIND(1D0)) :: irrfracbsoil
       REAL(KIND(1D0)) :: wetthresh
       TYPE(WATER_DIST_PRM) :: waterdist
    END TYPE LC_BSOIL_PRM
 
    TYPE, PUBLIC :: LC_WATER_PRM
       REAL(KIND(1D0)) :: sfr
       REAL(KIND(1D0)) :: emis
       TYPE(OHM_PRM) :: ohm
       TYPE(SOIL_PRM) :: soil
       REAL(KIND(1D0)) :: statelimit
       REAL(KIND(1D0)) :: irrfracwater
       REAL(KIND(1D0)) :: wetthresh ! ******* dummy variable *******
       REAL(KIND(1D0)) :: flowchange ! special term in water
    END TYPE LC_WATER_PRM
 
    ! ********** SUEWS_stateVars schema **********
    TYPE, PUBLIC :: anthroHEAT_STATE
       REAL(KIND(1D0)), DIMENSION(12) :: HDD_id !Heating Degree Days [degC d]
    END TYPE anthroHEAT_STATE
 
    TYPE, PUBLIC :: HYDRO_STATE
       ! REAL(KIND(1D0)) :: runofftowater   ! Fraction of above-ground runoff flowing to water surface during flooding
       REAL(KIND(1D0)), DIMENSION(nsurf) :: soilstore_surf ! Initial water stored in soil beneath `Bldgs` surface
       REAL(KIND(1D0)), DIMENSION(nsurf) :: state_surf ! Initial wetness condition on SUEWS land covers.
       REAL(KIND(1D0)), DIMENSION(9) :: WUDay_id ! Daily water use for EveTr, DecTr, Grass [mm]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: soilstore_roof ! Soil moisture of roof [mm]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: state_roof ! wetness status of roof [mm]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: soilstore_wall ! Soil moisture of wall [mm]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: state_wall ! wetness status of wall [mm]
    END TYPE HYDRO_STATE
 
    TYPE, PUBLIC :: HEAT_STATE
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: temp_roof ! interface temperature between depth layers in roof [degC]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: temp_wall ! interface temperature between depth layers in wall [degC]
       REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: temp_surf ! interface temperature between depth layers [degC]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: tsfc_roof ! roof surface temperature [degC]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: tsfc_wall ! wall surface temperature [degC]
       REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: tsfc_surf ! surface temperature [degC]
    END TYPE HEAT_STATE
 
    TYPE, PUBLIC :: OHM_STATE
       REAL(KIND(1D0)) :: qn_av ! weighted average of net all-wave radiation [W m-2]
       REAL(KIND(1D0)) :: dqndt ! rate of change of net radiation [W m-2 h-1]
       REAL(KIND(1D0)) :: qn_s_av ! weighted average of qn over snow [W m-2]
       REAL(KIND(1D0)) :: dqnsdt ! Rate of change of net radiation [W m-2 h-1]
    END TYPE OHM_STATE
 
    TYPE, PUBLIC :: PHENOLOGY_STATE
       REAL(KIND(1D0)), DIMENSION(NSURF) :: alb
       REAL(KIND(1D0)), DIMENSION(nvegsurf) :: lai_id ! Initial LAI values.
       REAL(KIND(1D0)), DIMENSION(nvegsurf) :: GDD_id ! Growing Degree Days [degC](see SUEWS_DailyState.f95)
       REAL(KIND(1D0)), DIMENSION(nvegsurf) :: SDD_id ! Senescence Degree Days [degC](see SUEWS_DailyState.f95)
       REAL(KIND(1D0)) :: porosity_id !
       REAL(KIND(1D0)) :: decidcap_id ! Storage capacity of deciduous surface `DecTr`; updated each day in simulaiton due to changes in LAI.
       REAL(KIND(1D0)) :: albDecTr_id !
       REAL(KIND(1D0)) :: albEveTr_id !
       REAL(KIND(1D0)) :: albGrass_id !
       REAL(KIND(1D0)) :: Tmin_id ! Daily minimum temperature [degC]
       REAL(KIND(1D0)) :: Tmax_id ! Daily maximum temperature [degC]
       REAL(KIND(1D0)) :: lenDay_id ! daytime length [h]
       REAL(KIND(1D0)), DIMENSION(6, NSURF) :: StoreDrainPrm ! coefficients used in drainage calculation [-]
    END TYPE PHENOLOGY_STATE
 
    TYPE, PUBLIC :: SNOW_STATE
       REAL(KIND(1D0)) :: snowfallCum !
       REAL(KIND(1D0)) :: snowalb ! albedo of snow [-]
       REAL(KIND(1D0)), DIMENSION(nsurf) :: icefrac ! fraction of ice in snowpack [-]
       REAL(KIND(1D0)), DIMENSION(nsurf) :: snowdens ! snow density [kg m-3]
       REAL(KIND(1D0)), DIMENSION(nsurf) :: snowfrac ! snow fraction [-]
       REAL(KIND(1D0)), DIMENSION(nsurf) :: snowpack ! snow water equivalent on each land cover [mm]
       REAL(KIND(1D0)), DIMENSION(nsurf) :: snowwater ! snow water[mm]
    END TYPE SNOW_STATE
 
    ! ********** SUEWS_forcing schema **********
    TYPE, PUBLIC :: SUEWS_FORCING
       REAL(KIND(1D0)) :: kdown !
       REAL(KIND(1D0)) :: ldown !
       REAL(KIND(1D0)) :: RH !
       REAL(KIND(1D0)) :: pres !
       REAL(KIND(1D0)) :: Tair !
       REAL(KIND(1D0)) :: U !
       REAL(KIND(1D0)) :: rain !
       REAL(KIND(1D0)) :: Wuh !  external water use
       REAL(KIND(1D0)) :: fcld !
       REAL(KIND(1D0)) :: LAI_obs !
       REAL(KIND(1D0)) :: snowfrac !
       REAL(KIND(1D0)) :: xsmd !
       REAL(KIND(1D0)) :: qf_obs !
       REAL(KIND(1D0)) :: qn1_obs !
       REAL(KIND(1D0)) :: qs_obs !
       REAL(KIND(1D0)) :: temp_c !
    END TYPE SUEWS_FORCING
 
    TYPE, PUBLIC :: SUEWS_TIMER
       INTEGER :: id !
       INTEGER :: imin !
       INTEGER :: isec !
       INTEGER :: it !
       INTEGER :: iy !
       INTEGER :: tstep !
       INTEGER :: tstep_prev !
       INTEGER :: dt_since_start !
    END TYPE SUEWS_TIMER
 
 END MODULE SUEWS_DEF_DTS