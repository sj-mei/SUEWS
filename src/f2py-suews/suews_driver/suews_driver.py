"""
Module suews_driver


Defined at suews_ctrl_driver.fpp lines 10-4460

"""
from __future__ import print_function, absolute_import, division
import _suews_driver
import f90wrap.runtime
import logging
import numpy

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("suews_driver.config")
class config(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=config)
    
    
    Defined at suews_ctrl_driver.fpp lines 50-54
    
    """
    def __init__(self, handle=None):
        """
        self = Config()
        
        
        Defined at suews_ctrl_driver.fpp lines 50-54
        
        
        Returns
        -------
        this : Config
        	Object to be constructed
        
        
        Automatically generated constructor for config
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _suews_driver.f90wrap_config_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Config
        
        
        Defined at suews_ctrl_driver.fpp lines 50-54
        
        Parameters
        ----------
        this : Config
        	Object to be destructed
        
        
        Automatically generated destructor for config
        """
        if self._alloc:
            _suews_driver.f90wrap_config_finalise(this=self._handle)
    
    @property
    def var1(self):
        """
        Element var1 ftype=integer  pytype=int
        
        
        Defined at suews_ctrl_driver.fpp line 51
        
        """
        return _suews_driver.f90wrap_config__get__var1(self._handle)
    
    @var1.setter
    def var1(self, var1):
        _suews_driver.f90wrap_config__set__var1(self._handle, var1)
    
    @property
    def var2(self):
        """
        Element var2 ftype=integer  pytype=int
        
        
        Defined at suews_ctrl_driver.fpp line 52
        
        """
        return _suews_driver.f90wrap_config__get__var2(self._handle)
    
    @var2.setter
    def var2(self, var2):
        _suews_driver.f90wrap_config__set__var2(self._handle, var2)
    
    def __str__(self):
        ret = ['<config>{\n']
        ret.append('    var1 : ')
        ret.append(repr(self.var1))
        ret.append(',\n    var2 : ')
        ret.append(repr(self.var2))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def var2add_two(self):
    """
    res = var2add_two(self)
    
    
    Defined at suews_ctrl_driver.fpp lines 57-60
    
    Parameters
    ----------
    arg_type : Config
    
    Returns
    -------
    res : int
    
    """
    res = _suews_driver.f90wrap_var2add_two(arg_type=self._handle)
    return res

def suews_cal_main(ah_min, ahprof_24hr, ah_slope_cooling, ah_slope_heating, alb, \
    albmax_dectr, albmax_evetr, albmax_grass, albmin_dectr, albmin_evetr, \
    albmin_grass, alpha_bioco2, alpha_enh_bioco2, alt, kdown, avrh, avu1, baset, \
    basete, beta_bioco2, beta_enh_bioco2, bldgh, capmax_dec, capmin_dec, \
    chanohm, co2pointsource, cpanohm, crwmax, crwmin, daywat, daywatper, \
    dectreeh, diagmethod, diagnose, drainrt, dt_since_start, dqndt, qn_av, \
    dqnsdt, qn_s_av, ef_umolco2perj, emis, emissionsmethod, enef_v_jkm, enddls, \
    evetreeh, faibldg, faidectree, faievetree, faut, fcef_v_kgkm, fcld_obs, \
    flowchange, frfossilfuel_heat, frfossilfuel_nonheat, g_max, g_k, g_q_base, \
    g_q_shape, g_t, g_sm, gdd_id, gddfull, gridiv, gsmodel, h_maintain, hdd_id, \
    humactivity_24hr, icefrac, id, ie_a, ie_end, ie_m, ie_start, imin, \
    internalwateruse_h, irrfracpaved, irrfracbldgs, irrfracevetr, irrfracdectr, \
    irrfracgrass, irrfracbsoil, irrfracwater, isec, it, iy, kkanohm, kmax, \
    lai_id, laimax, laimin, lai_obs, laipower, laitype, lat, lenday_id, \
    ldown_obs, lng, maxconductance, maxfcmetab, maxqfmetab, snowwater, \
    metforcingdata_grid, minfcmetab, minqfmetab, min_res_bioco2, narp_emis_snow, \
    narp_trans_site, netradiationmethod, nlayer, n_vegetation_region_urban, \
    n_stream_sw_urban, n_stream_lw_urban, sw_dn_direct_frac, air_ext_sw, \
    air_ssa_sw, veg_ssa_sw, air_ext_lw, air_ssa_lw, veg_ssa_lw, veg_fsd_const, \
    veg_contact_fraction_const, ground_albedo_dir_mult_fact, \
    use_sw_direct_albedo, height, building_frac, veg_frac, building_scale, \
    veg_scale, alb_roof, emis_roof, alb_wall, emis_wall, \
    roof_albedo_dir_mult_fact, wall_specular_frac, ohm_coef, ohmincqf, \
    ohm_threshsw, ohm_threshwd, pipecapacity, popdensdaytime, popdensnighttime, \
    popprof_24hr, pormax_dec, pormin_dec, precip, preciplimit, preciplimitalb, \
    press_hpa, qf0_beu, qf_a, qf_b, qf_c, qn1_obs, qs_obs, qf_obs, radmeltfact, \
    raincover, rainmaxres, resp_a, resp_b, roughlenheatmethod, \
    roughlenmommethod, runofftowater, s1, s2, sathydraulicconduct, sddfull, \
    sdd_id, smdmethod, snowalb, snowalbmax, snowalbmin, snowpacklimit, snowdens, \
    snowdensmax, snowdensmin, snowfallcum, snowfrac, snowlimbldg, snowlimpaved, \
    snowfrac_obs, snowpack, snowprof_24hr, snowuse, soildepth, stabilitymethod, \
    startdls, soilstore_surf, soilstorecap_surf, state_surf, statelimit_surf, \
    wetthresh_surf, soilstore_roof, soilstorecap_roof, state_roof, \
    statelimit_roof, wetthresh_roof, soilstore_wall, soilstorecap_wall, \
    state_wall, statelimit_wall, wetthresh_wall, storageheatmethod, \
    storedrainprm, surfacearea, tair_av, tau_a, tau_f, tau_r, tmax_id, tmin_id, \
    baset_cooling, baset_heating, temp_c, tempmeltfact, th, theta_bioco2, \
    timezone, tl, trafficrate, trafficunits, sfr_surf, tsfc_roof, tsfc_wall, \
    tsfc_surf, temp_roof, temp_wall, temp_surf, tin_roof, tin_wall, tin_surf, \
    k_roof, k_wall, k_surf, cp_roof, cp_wall, cp_surf, dz_roof, dz_wall, \
    dz_surf, traffprof_24hr, ts5mindata_ir, tstep, tstep_prev, veg_type, \
    waterdist, waterusemethod, wu_m3, wuday_id, decidcap_id, albdectr_id, \
    albevetr_id, albgrass_id, porosity_id, wuprofa_24hr, wuprofm_24hr, xsmd, z, \
    z0m_in, zdm_in, datetimeline, dataoutlinesuews, dataoutlinesnow, \
    dataoutlineestm, dataoutlinersl, dataoutlinebeers, dataoutlinedebug, \
    dataoutlinespartacus, dataoutlineestmext, dailystateline):
    """
    suews_cal_main(ah_min, ahprof_24hr, ah_slope_cooling, ah_slope_heating, alb, \
        albmax_dectr, albmax_evetr, albmax_grass, albmin_dectr, albmin_evetr, \
        albmin_grass, alpha_bioco2, alpha_enh_bioco2, alt, kdown, avrh, avu1, baset, \
        basete, beta_bioco2, beta_enh_bioco2, bldgh, capmax_dec, capmin_dec, \
        chanohm, co2pointsource, cpanohm, crwmax, crwmin, daywat, daywatper, \
        dectreeh, diagmethod, diagnose, drainrt, dt_since_start, dqndt, qn_av, \
        dqnsdt, qn_s_av, ef_umolco2perj, emis, emissionsmethod, enef_v_jkm, enddls, \
        evetreeh, faibldg, faidectree, faievetree, faut, fcef_v_kgkm, fcld_obs, \
        flowchange, frfossilfuel_heat, frfossilfuel_nonheat, g_max, g_k, g_q_base, \
        g_q_shape, g_t, g_sm, gdd_id, gddfull, gridiv, gsmodel, h_maintain, hdd_id, \
        humactivity_24hr, icefrac, id, ie_a, ie_end, ie_m, ie_start, imin, \
        internalwateruse_h, irrfracpaved, irrfracbldgs, irrfracevetr, irrfracdectr, \
        irrfracgrass, irrfracbsoil, irrfracwater, isec, it, iy, kkanohm, kmax, \
        lai_id, laimax, laimin, lai_obs, laipower, laitype, lat, lenday_id, \
        ldown_obs, lng, maxconductance, maxfcmetab, maxqfmetab, snowwater, \
        metforcingdata_grid, minfcmetab, minqfmetab, min_res_bioco2, narp_emis_snow, \
        narp_trans_site, netradiationmethod, nlayer, n_vegetation_region_urban, \
        n_stream_sw_urban, n_stream_lw_urban, sw_dn_direct_frac, air_ext_sw, \
        air_ssa_sw, veg_ssa_sw, air_ext_lw, air_ssa_lw, veg_ssa_lw, veg_fsd_const, \
        veg_contact_fraction_const, ground_albedo_dir_mult_fact, \
        use_sw_direct_albedo, height, building_frac, veg_frac, building_scale, \
        veg_scale, alb_roof, emis_roof, alb_wall, emis_wall, \
        roof_albedo_dir_mult_fact, wall_specular_frac, ohm_coef, ohmincqf, \
        ohm_threshsw, ohm_threshwd, pipecapacity, popdensdaytime, popdensnighttime, \
        popprof_24hr, pormax_dec, pormin_dec, precip, preciplimit, preciplimitalb, \
        press_hpa, qf0_beu, qf_a, qf_b, qf_c, qn1_obs, qs_obs, qf_obs, radmeltfact, \
        raincover, rainmaxres, resp_a, resp_b, roughlenheatmethod, \
        roughlenmommethod, runofftowater, s1, s2, sathydraulicconduct, sddfull, \
        sdd_id, smdmethod, snowalb, snowalbmax, snowalbmin, snowpacklimit, snowdens, \
        snowdensmax, snowdensmin, snowfallcum, snowfrac, snowlimbldg, snowlimpaved, \
        snowfrac_obs, snowpack, snowprof_24hr, snowuse, soildepth, stabilitymethod, \
        startdls, soilstore_surf, soilstorecap_surf, state_surf, statelimit_surf, \
        wetthresh_surf, soilstore_roof, soilstorecap_roof, state_roof, \
        statelimit_roof, wetthresh_roof, soilstore_wall, soilstorecap_wall, \
        state_wall, statelimit_wall, wetthresh_wall, storageheatmethod, \
        storedrainprm, surfacearea, tair_av, tau_a, tau_f, tau_r, tmax_id, tmin_id, \
        baset_cooling, baset_heating, temp_c, tempmeltfact, th, theta_bioco2, \
        timezone, tl, trafficrate, trafficunits, sfr_surf, tsfc_roof, tsfc_wall, \
        tsfc_surf, temp_roof, temp_wall, temp_surf, tin_roof, tin_wall, tin_surf, \
        k_roof, k_wall, k_surf, cp_roof, cp_wall, cp_surf, dz_roof, dz_wall, \
        dz_surf, traffprof_24hr, ts5mindata_ir, tstep, tstep_prev, veg_type, \
        waterdist, waterusemethod, wu_m3, wuday_id, decidcap_id, albdectr_id, \
        albevetr_id, albgrass_id, porosity_id, wuprofa_24hr, wuprofm_24hr, xsmd, z, \
        z0m_in, zdm_in, datetimeline, dataoutlinesuews, dataoutlinesnow, \
        dataoutlineestm, dataoutlinersl, dataoutlinebeers, dataoutlinedebug, \
        dataoutlinespartacus, dataoutlineestmext, dailystateline)
    
    
    Defined at suews_ctrl_driver.fpp lines 131-1438
    
    Parameters
    ----------
    ah_min : float array
    ahprof_24hr : float array
    ah_slope_cooling : float array
    ah_slope_heating : float array
    alb : float array
    albmax_dectr : float
    albmax_evetr : float
    albmax_grass : float
    albmin_dectr : float
    albmin_evetr : float
    albmin_grass : float
    alpha_bioco2 : float array
    alpha_enh_bioco2 : float array
    alt : float
    kdown : float
    avrh : float
    avu1 : float
    baset : float array
    basete : float array
    beta_bioco2 : float array
    beta_enh_bioco2 : float array
    bldgh : float
    capmax_dec : float
    capmin_dec : float
    chanohm : float array
    co2pointsource : float
    cpanohm : float array
    crwmax : float
    crwmin : float
    daywat : float array
    daywatper : float array
    dectreeh : float
    diagmethod : int
    diagnose : int
    drainrt : float
    dt_since_start : int
    dqndt : float
    qn_av : float
    dqnsdt : float
    qn_s_av : float
    ef_umolco2perj : float
    emis : float array
    emissionsmethod : int
    enef_v_jkm : float
    enddls : int
    evetreeh : float
    faibldg : float
    faidectree : float
    faievetree : float
    faut : float
    fcef_v_kgkm : float array
    fcld_obs : float
    flowchange : float
    frfossilfuel_heat : float
    frfossilfuel_nonheat : float
    g_max : float
    g_k : float
    g_q_base : float
    g_q_shape : float
    g_t : float
    g_sm : float
    gdd_id : float array
    gddfull : float array
    gridiv : int
    gsmodel : int
    h_maintain : float
    hdd_id : float array
    humactivity_24hr : float array
    icefrac : float array
    id : int
    ie_a : float array
    ie_end : int
    ie_m : float array
    ie_start : int
    imin : int
    internalwateruse_h : float
    irrfracpaved : float
    irrfracbldgs : float
    irrfracevetr : float
    irrfracdectr : float
    irrfracgrass : float
    irrfracbsoil : float
    irrfracwater : float
    isec : int
    it : int
    iy : int
    kkanohm : float array
    kmax : float
    lai_id : float array
    laimax : float array
    laimin : float array
    lai_obs : float
    laipower : float array
    laitype : int array
    lat : float
    lenday_id : float
    ldown_obs : float
    lng : float
    maxconductance : float array
    maxfcmetab : float
    maxqfmetab : float
    snowwater : float array
    metforcingdata_grid : float array
    minfcmetab : float
    minqfmetab : float
    min_res_bioco2 : float array
    narp_emis_snow : float
    narp_trans_site : float
    netradiationmethod : int
    nlayer : int
    n_vegetation_region_urban : int
    n_stream_sw_urban : int
    n_stream_lw_urban : int
    sw_dn_direct_frac : float
    air_ext_sw : float
    air_ssa_sw : float
    veg_ssa_sw : float
    air_ext_lw : float
    air_ssa_lw : float
    veg_ssa_lw : float
    veg_fsd_const : float
    veg_contact_fraction_const : float
    ground_albedo_dir_mult_fact : float
    use_sw_direct_albedo : bool
    height : float array
    building_frac : float array
    veg_frac : float array
    building_scale : float array
    veg_scale : float array
    alb_roof : float array
    emis_roof : float array
    alb_wall : float array
    emis_wall : float array
    roof_albedo_dir_mult_fact : float array
    wall_specular_frac : float array
    ohm_coef : float array
    ohmincqf : int
    ohm_threshsw : float array
    ohm_threshwd : float array
    pipecapacity : float
    popdensdaytime : float array
    popdensnighttime : float
    popprof_24hr : float array
    pormax_dec : float
    pormin_dec : float
    precip : float
    preciplimit : float
    preciplimitalb : float
    press_hpa : float
    qf0_beu : float array
    qf_a : float array
    qf_b : float array
    qf_c : float array
    qn1_obs : float
    qs_obs : float
    qf_obs : float
    radmeltfact : float
    raincover : float
    rainmaxres : float
    resp_a : float array
    resp_b : float array
    roughlenheatmethod : int
    roughlenmommethod : int
    runofftowater : float
    s1 : float
    s2 : float
    sathydraulicconduct : float array
    sddfull : float array
    sdd_id : float array
    smdmethod : int
    snowalb : float
    snowalbmax : float
    snowalbmin : float
    snowpacklimit : float array
    snowdens : float array
    snowdensmax : float
    snowdensmin : float
    snowfallcum : float
    snowfrac : float array
    snowlimbldg : float
    snowlimpaved : float
    snowfrac_obs : float
    snowpack : float array
    snowprof_24hr : float array
    snowuse : int
    soildepth : float array
    stabilitymethod : int
    startdls : int
    soilstore_surf : float array
    soilstorecap_surf : float array
    state_surf : float array
    statelimit_surf : float array
    wetthresh_surf : float array
    soilstore_roof : float array
    soilstorecap_roof : float array
    state_roof : float array
    statelimit_roof : float array
    wetthresh_roof : float array
    soilstore_wall : float array
    soilstorecap_wall : float array
    state_wall : float array
    statelimit_wall : float array
    wetthresh_wall : float array
    storageheatmethod : int
    storedrainprm : float array
    surfacearea : float
    tair_av : float
    tau_a : float
    tau_f : float
    tau_r : float
    tmax_id : float
    tmin_id : float
    baset_cooling : float array
    baset_heating : float array
    temp_c : float
    tempmeltfact : float
    th : float
    theta_bioco2 : float array
    timezone : float
    tl : float
    trafficrate : float array
    trafficunits : float
    sfr_surf : float array
    tsfc_roof : float array
    tsfc_wall : float array
    tsfc_surf : float array
    temp_roof : float array
    temp_wall : float array
    temp_surf : float array
    tin_roof : float array
    tin_wall : float array
    tin_surf : float array
    k_roof : float array
    k_wall : float array
    k_surf : float array
    cp_roof : float array
    cp_wall : float array
    cp_surf : float array
    dz_roof : float array
    dz_wall : float array
    dz_surf : float array
    traffprof_24hr : float array
    ts5mindata_ir : float array
    tstep : int
    tstep_prev : int
    veg_type : int
    waterdist : float array
    waterusemethod : int
    wu_m3 : float
    wuday_id : float array
    decidcap_id : float
    albdectr_id : float
    albevetr_id : float
    albgrass_id : float
    porosity_id : float
    wuprofa_24hr : float array
    wuprofm_24hr : float array
    xsmd : float
    z : float
    z0m_in : float
    zdm_in : float
    datetimeline : float array
    dataoutlinesuews : float array
    dataoutlinesnow : float array
    dataoutlineestm : float array
    dataoutlinersl : float array
    dataoutlinebeers : float array
    dataoutlinedebug : float array
    dataoutlinespartacus : float array
    dataoutlineestmext : float array
    dailystateline : float array
    
    ==============main calculation start=======================
    ==============surface roughness calculation=======================
    """
    _suews_driver.f90wrap_suews_cal_main(ah_min=ah_min, ahprof_24hr=ahprof_24hr, \
        ah_slope_cooling=ah_slope_cooling, ah_slope_heating=ah_slope_heating, \
        alb=alb, albmax_dectr=albmax_dectr, albmax_evetr=albmax_evetr, \
        albmax_grass=albmax_grass, albmin_dectr=albmin_dectr, \
        albmin_evetr=albmin_evetr, albmin_grass=albmin_grass, \
        alpha_bioco2=alpha_bioco2, alpha_enh_bioco2=alpha_enh_bioco2, alt=alt, \
        kdown=kdown, avrh=avrh, avu1=avu1, baset=baset, basete=basete, \
        beta_bioco2=beta_bioco2, beta_enh_bioco2=beta_enh_bioco2, bldgh=bldgh, \
        capmax_dec=capmax_dec, capmin_dec=capmin_dec, chanohm=chanohm, \
        co2pointsource=co2pointsource, cpanohm=cpanohm, crwmax=crwmax, \
        crwmin=crwmin, daywat=daywat, daywatper=daywatper, dectreeh=dectreeh, \
        diagmethod=diagmethod, diagnose=diagnose, drainrt=drainrt, \
        dt_since_start=dt_since_start, dqndt=dqndt, qn_av=qn_av, dqnsdt=dqnsdt, \
        qn_s_av=qn_s_av, ef_umolco2perj=ef_umolco2perj, emis=emis, \
        emissionsmethod=emissionsmethod, enef_v_jkm=enef_v_jkm, enddls=enddls, \
        evetreeh=evetreeh, faibldg=faibldg, faidectree=faidectree, \
        faievetree=faievetree, faut=faut, fcef_v_kgkm=fcef_v_kgkm, \
        fcld_obs=fcld_obs, flowchange=flowchange, \
        frfossilfuel_heat=frfossilfuel_heat, \
        frfossilfuel_nonheat=frfossilfuel_nonheat, g_max=g_max, g_k=g_k, \
        g_q_base=g_q_base, g_q_shape=g_q_shape, g_t=g_t, g_sm=g_sm, gdd_id=gdd_id, \
        gddfull=gddfull, gridiv=gridiv, gsmodel=gsmodel, h_maintain=h_maintain, \
        hdd_id=hdd_id, humactivity_24hr=humactivity_24hr, icefrac=icefrac, id=id, \
        ie_a=ie_a, ie_end=ie_end, ie_m=ie_m, ie_start=ie_start, imin=imin, \
        internalwateruse_h=internalwateruse_h, irrfracpaved=irrfracpaved, \
        irrfracbldgs=irrfracbldgs, irrfracevetr=irrfracevetr, \
        irrfracdectr=irrfracdectr, irrfracgrass=irrfracgrass, \
        irrfracbsoil=irrfracbsoil, irrfracwater=irrfracwater, isec=isec, it=it, \
        iy=iy, kkanohm=kkanohm, kmax=kmax, lai_id=lai_id, laimax=laimax, \
        laimin=laimin, lai_obs=lai_obs, laipower=laipower, laitype=laitype, lat=lat, \
        lenday_id=lenday_id, ldown_obs=ldown_obs, lng=lng, \
        maxconductance=maxconductance, maxfcmetab=maxfcmetab, maxqfmetab=maxqfmetab, \
        snowwater=snowwater, metforcingdata_grid=metforcingdata_grid, \
        minfcmetab=minfcmetab, minqfmetab=minqfmetab, min_res_bioco2=min_res_bioco2, \
        narp_emis_snow=narp_emis_snow, narp_trans_site=narp_trans_site, \
        netradiationmethod=netradiationmethod, nlayer=nlayer, \
        n_vegetation_region_urban=n_vegetation_region_urban, \
        n_stream_sw_urban=n_stream_sw_urban, n_stream_lw_urban=n_stream_lw_urban, \
        sw_dn_direct_frac=sw_dn_direct_frac, air_ext_sw=air_ext_sw, \
        air_ssa_sw=air_ssa_sw, veg_ssa_sw=veg_ssa_sw, air_ext_lw=air_ext_lw, \
        air_ssa_lw=air_ssa_lw, veg_ssa_lw=veg_ssa_lw, veg_fsd_const=veg_fsd_const, \
        veg_contact_fraction_const=veg_contact_fraction_const, \
        ground_albedo_dir_mult_fact=ground_albedo_dir_mult_fact, \
        use_sw_direct_albedo=use_sw_direct_albedo, height=height, \
        building_frac=building_frac, veg_frac=veg_frac, \
        building_scale=building_scale, veg_scale=veg_scale, alb_roof=alb_roof, \
        emis_roof=emis_roof, alb_wall=alb_wall, emis_wall=emis_wall, \
        roof_albedo_dir_mult_fact=roof_albedo_dir_mult_fact, \
        wall_specular_frac=wall_specular_frac, ohm_coef=ohm_coef, ohmincqf=ohmincqf, \
        ohm_threshsw=ohm_threshsw, ohm_threshwd=ohm_threshwd, \
        pipecapacity=pipecapacity, popdensdaytime=popdensdaytime, \
        popdensnighttime=popdensnighttime, popprof_24hr=popprof_24hr, \
        pormax_dec=pormax_dec, pormin_dec=pormin_dec, precip=precip, \
        preciplimit=preciplimit, preciplimitalb=preciplimitalb, press_hpa=press_hpa, \
        qf0_beu=qf0_beu, qf_a=qf_a, qf_b=qf_b, qf_c=qf_c, qn1_obs=qn1_obs, \
        qs_obs=qs_obs, qf_obs=qf_obs, radmeltfact=radmeltfact, raincover=raincover, \
        rainmaxres=rainmaxres, resp_a=resp_a, resp_b=resp_b, \
        roughlenheatmethod=roughlenheatmethod, roughlenmommethod=roughlenmommethod, \
        runofftowater=runofftowater, s1=s1, s2=s2, \
        sathydraulicconduct=sathydraulicconduct, sddfull=sddfull, sdd_id=sdd_id, \
        smdmethod=smdmethod, snowalb=snowalb, snowalbmax=snowalbmax, \
        snowalbmin=snowalbmin, snowpacklimit=snowpacklimit, snowdens=snowdens, \
        snowdensmax=snowdensmax, snowdensmin=snowdensmin, snowfallcum=snowfallcum, \
        snowfrac=snowfrac, snowlimbldg=snowlimbldg, snowlimpaved=snowlimpaved, \
        snowfrac_obs=snowfrac_obs, snowpack=snowpack, snowprof_24hr=snowprof_24hr, \
        snowuse=snowuse, soildepth=soildepth, stabilitymethod=stabilitymethod, \
        startdls=startdls, soilstore_surf=soilstore_surf, \
        soilstorecap_surf=soilstorecap_surf, state_surf=state_surf, \
        statelimit_surf=statelimit_surf, wetthresh_surf=wetthresh_surf, \
        soilstore_roof=soilstore_roof, soilstorecap_roof=soilstorecap_roof, \
        state_roof=state_roof, statelimit_roof=statelimit_roof, \
        wetthresh_roof=wetthresh_roof, soilstore_wall=soilstore_wall, \
        soilstorecap_wall=soilstorecap_wall, state_wall=state_wall, \
        statelimit_wall=statelimit_wall, wetthresh_wall=wetthresh_wall, \
        storageheatmethod=storageheatmethod, storedrainprm=storedrainprm, \
        surfacearea=surfacearea, tair_av=tair_av, tau_a=tau_a, tau_f=tau_f, \
        tau_r=tau_r, tmax_id=tmax_id, tmin_id=tmin_id, baset_cooling=baset_cooling, \
        baset_heating=baset_heating, temp_c=temp_c, tempmeltfact=tempmeltfact, \
        th=th, theta_bioco2=theta_bioco2, timezone=timezone, tl=tl, \
        trafficrate=trafficrate, trafficunits=trafficunits, sfr_surf=sfr_surf, \
        tsfc_roof=tsfc_roof, tsfc_wall=tsfc_wall, tsfc_surf=tsfc_surf, \
        temp_roof=temp_roof, temp_wall=temp_wall, temp_surf=temp_surf, \
        tin_roof=tin_roof, tin_wall=tin_wall, tin_surf=tin_surf, k_roof=k_roof, \
        k_wall=k_wall, k_surf=k_surf, cp_roof=cp_roof, cp_wall=cp_wall, \
        cp_surf=cp_surf, dz_roof=dz_roof, dz_wall=dz_wall, dz_surf=dz_surf, \
        traffprof_24hr=traffprof_24hr, ts5mindata_ir=ts5mindata_ir, tstep=tstep, \
        tstep_prev=tstep_prev, veg_type=veg_type, waterdist=waterdist, \
        waterusemethod=waterusemethod, wu_m3=wu_m3, wuday_id=wuday_id, \
        decidcap_id=decidcap_id, albdectr_id=albdectr_id, albevetr_id=albevetr_id, \
        albgrass_id=albgrass_id, porosity_id=porosity_id, wuprofa_24hr=wuprofa_24hr, \
        wuprofm_24hr=wuprofm_24hr, xsmd=xsmd, z=z, z0m_in=z0m_in, zdm_in=zdm_in, \
        datetimeline=datetimeline, dataoutlinesuews=dataoutlinesuews, \
        dataoutlinesnow=dataoutlinesnow, dataoutlineestm=dataoutlineestm, \
        dataoutlinersl=dataoutlinersl, dataoutlinebeers=dataoutlinebeers, \
        dataoutlinedebug=dataoutlinedebug, \
        dataoutlinespartacus=dataoutlinespartacus, \
        dataoutlineestmext=dataoutlineestmext, dailystateline=dailystateline)

def suews_cal_anthropogenicemission(ah_min, ahprof_24hr, ah_slope_cooling, \
    ah_slope_heating, co2pointsource, dayofweek_id, dls, ef_umolco2perj, \
    emissionsmethod, enef_v_jkm, fcef_v_kgkm, frfossilfuel_heat, \
    frfossilfuel_nonheat, hdd_id, humactivity_24hr, imin, it, maxfcmetab, \
    maxqfmetab, minfcmetab, minqfmetab, popdensdaytime, popdensnighttime, \
    popprof_24hr, qf0_beu, qf_a, qf_b, qf_c, qf_obs, surfacearea, baset_cooling, \
    baset_heating, temp_c, trafficrate, trafficunits, traffprof_24hr):
    """
    qf, qf_sahp, fc_anthro, fc_build, fc_metab, fc_point, fc_traff = \
        suews_cal_anthropogenicemission(ah_min, ahprof_24hr, ah_slope_cooling, \
        ah_slope_heating, co2pointsource, dayofweek_id, dls, ef_umolco2perj, \
        emissionsmethod, enef_v_jkm, fcef_v_kgkm, frfossilfuel_heat, \
        frfossilfuel_nonheat, hdd_id, humactivity_24hr, imin, it, maxfcmetab, \
        maxqfmetab, minfcmetab, minqfmetab, popdensdaytime, popdensnighttime, \
        popprof_24hr, qf0_beu, qf_a, qf_b, qf_c, qf_obs, surfacearea, baset_cooling, \
        baset_heating, temp_c, trafficrate, trafficunits, traffprof_24hr)
    
    
    Defined at suews_ctrl_driver.fpp lines 1450-1530
    
    Parameters
    ----------
    ah_min : float array
    ahprof_24hr : float array
    ah_slope_cooling : float array
    ah_slope_heating : float array
    co2pointsource : float
    dayofweek_id : int array
    dls : int
    ef_umolco2perj : float
    emissionsmethod : int
    enef_v_jkm : float
    fcef_v_kgkm : float array
    frfossilfuel_heat : float
    frfossilfuel_nonheat : float
    hdd_id : float array
    humactivity_24hr : float array
    imin : int
    it : int
    maxfcmetab : float
    maxqfmetab : float
    minfcmetab : float
    minqfmetab : float
    popdensdaytime : float array
    popdensnighttime : float
    popprof_24hr : float array
    qf0_beu : float array
    qf_a : float array
    qf_b : float array
    qf_c : float array
    qf_obs : float
    surfacearea : float
    baset_cooling : float array
    baset_heating : float array
    temp_c : float
    trafficrate : float array
    trafficunits : float
    traffprof_24hr : float array
    
    Returns
    -------
    qf : float
    qf_sahp : float
    fc_anthro : float
    fc_build : float
    fc_metab : float
    fc_point : float
    fc_traff : float
    
    """
    qf, qf_sahp, fc_anthro, fc_build, fc_metab, fc_point, fc_traff = \
        _suews_driver.f90wrap_suews_cal_anthropogenicemission(ah_min=ah_min, \
        ahprof_24hr=ahprof_24hr, ah_slope_cooling=ah_slope_cooling, \
        ah_slope_heating=ah_slope_heating, co2pointsource=co2pointsource, \
        dayofweek_id=dayofweek_id, dls=dls, ef_umolco2perj=ef_umolco2perj, \
        emissionsmethod=emissionsmethod, enef_v_jkm=enef_v_jkm, \
        fcef_v_kgkm=fcef_v_kgkm, frfossilfuel_heat=frfossilfuel_heat, \
        frfossilfuel_nonheat=frfossilfuel_nonheat, hdd_id=hdd_id, \
        humactivity_24hr=humactivity_24hr, imin=imin, it=it, maxfcmetab=maxfcmetab, \
        maxqfmetab=maxqfmetab, minfcmetab=minfcmetab, minqfmetab=minqfmetab, \
        popdensdaytime=popdensdaytime, popdensnighttime=popdensnighttime, \
        popprof_24hr=popprof_24hr, qf0_beu=qf0_beu, qf_a=qf_a, qf_b=qf_b, qf_c=qf_c, \
        qf_obs=qf_obs, surfacearea=surfacearea, baset_cooling=baset_cooling, \
        baset_heating=baset_heating, temp_c=temp_c, trafficrate=trafficrate, \
        trafficunits=trafficunits, traffprof_24hr=traffprof_24hr)
    return qf, qf_sahp, fc_anthro, fc_build, fc_metab, fc_point, fc_traff

def suews_cal_biogenco2(alpha_bioco2, alpha_enh_bioco2, avkdn, avrh, \
    beta_bioco2, beta_enh_bioco2, dectime, diagnose, emissionsmethod, fc_anthro, \
    g_max, g_k, g_q_base, g_q_shape, g_t, g_sm, gfunc, gsmodel, id, it, kmax, \
    lai_id, laimin, laimax, maxconductance, min_res_bioco2, press_hpa, resp_a, \
    resp_b, s1, s2, sfr_surf, smdmethod, snowfrac, t2_c, temp_c, theta_bioco2, \
    th, tl, vsmd, xsmd):
    """
    fc, fc_biogen, fc_photo, fc_respi = suews_cal_biogenco2(alpha_bioco2, \
        alpha_enh_bioco2, avkdn, avrh, beta_bioco2, beta_enh_bioco2, dectime, \
        diagnose, emissionsmethod, fc_anthro, g_max, g_k, g_q_base, g_q_shape, g_t, \
        g_sm, gfunc, gsmodel, id, it, kmax, lai_id, laimin, laimax, maxconductance, \
        min_res_bioco2, press_hpa, resp_a, resp_b, s1, s2, sfr_surf, smdmethod, \
        snowfrac, t2_c, temp_c, theta_bioco2, th, tl, vsmd, xsmd)
    
    
    Defined at suews_ctrl_driver.fpp lines 1540-1645
    
    Parameters
    ----------
    alpha_bioco2 : float array
    alpha_enh_bioco2 : float array
    avkdn : float
    avrh : float
    beta_bioco2 : float array
    beta_enh_bioco2 : float array
    dectime : float
    diagnose : int
    emissionsmethod : int
    fc_anthro : float
    g_max : float
    g_k : float
    g_q_base : float
    g_q_shape : float
    g_t : float
    g_sm : float
    gfunc : float
    gsmodel : int
    id : int
    it : int
    kmax : float
    lai_id : float array
    laimin : float array
    laimax : float array
    maxconductance : float array
    min_res_bioco2 : float array
    press_hpa : float
    resp_a : float array
    resp_b : float array
    s1 : float
    s2 : float
    sfr_surf : float array
    smdmethod : int
    snowfrac : float array
    t2_c : float
    temp_c : float
    theta_bioco2 : float array
    th : float
    tl : float
    vsmd : float
    xsmd : float
    
    Returns
    -------
    fc : float
    fc_biogen : float
    fc_photo : float
    fc_respi : float
    
    """
    fc, fc_biogen, fc_photo, fc_respi = \
        _suews_driver.f90wrap_suews_cal_biogenco2(alpha_bioco2=alpha_bioco2, \
        alpha_enh_bioco2=alpha_enh_bioco2, avkdn=avkdn, avrh=avrh, \
        beta_bioco2=beta_bioco2, beta_enh_bioco2=beta_enh_bioco2, dectime=dectime, \
        diagnose=diagnose, emissionsmethod=emissionsmethod, fc_anthro=fc_anthro, \
        g_max=g_max, g_k=g_k, g_q_base=g_q_base, g_q_shape=g_q_shape, g_t=g_t, \
        g_sm=g_sm, gfunc=gfunc, gsmodel=gsmodel, id=id, it=it, kmax=kmax, \
        lai_id=lai_id, laimin=laimin, laimax=laimax, maxconductance=maxconductance, \
        min_res_bioco2=min_res_bioco2, press_hpa=press_hpa, resp_a=resp_a, \
        resp_b=resp_b, s1=s1, s2=s2, sfr_surf=sfr_surf, smdmethod=smdmethod, \
        snowfrac=snowfrac, t2_c=t2_c, temp_c=temp_c, theta_bioco2=theta_bioco2, \
        th=th, tl=tl, vsmd=vsmd, xsmd=xsmd)
    return fc, fc_biogen, fc_photo, fc_respi

def suews_cal_qn(storageheatmethod, netradiationmethod, snowuse, tstep, nlayer, \
    snowpack_prev, tau_a, tau_f, snowalbmax, snowalbmin, diagnose, ldown_obs, \
    fcld_obs, dectime, zenith_deg, tsurf_0, kdown, tair_c, avrh, ea_hpa, \
    qn1_obs, snowalb_prev, snowfrac_prev, diagqn, narp_trans_site, \
    narp_emis_snow, icefrac, sfr_surf, sfr_roof, sfr_wall, tsfc_surf, tsfc_roof, \
    tsfc_wall, emis, alb_prev, albdectr_id, albevetr_id, albgrass_id, lai_id, \
    n_vegetation_region_urban, n_stream_sw_urban, n_stream_lw_urban, \
    sw_dn_direct_frac, air_ext_sw, air_ssa_sw, veg_ssa_sw, air_ext_lw, \
    air_ssa_lw, veg_ssa_lw, veg_fsd_const, veg_contact_fraction_const, \
    ground_albedo_dir_mult_fact, use_sw_direct_albedo, height, building_frac, \
    veg_frac, building_scale, veg_scale, alb_roof, emis_roof, alb_wall, \
    emis_wall, roof_albedo_dir_mult_fact, wall_specular_frac, alb_next, qn_surf, \
    qn_roof, qn_wall, qn_ind_snow, kup_ind_snow, tsurf_ind_snow, tsurf_ind, \
    dataoutlinespartacus):
    """
    ldown, fcld, qn, qn_snowfree, qn_snow, kclear, kup, lup, tsurf, albedo_snow, \
        snowalb_next = suews_cal_qn(storageheatmethod, netradiationmethod, snowuse, \
        tstep, nlayer, snowpack_prev, tau_a, tau_f, snowalbmax, snowalbmin, \
        diagnose, ldown_obs, fcld_obs, dectime, zenith_deg, tsurf_0, kdown, tair_c, \
        avrh, ea_hpa, qn1_obs, snowalb_prev, snowfrac_prev, diagqn, narp_trans_site, \
        narp_emis_snow, icefrac, sfr_surf, sfr_roof, sfr_wall, tsfc_surf, tsfc_roof, \
        tsfc_wall, emis, alb_prev, albdectr_id, albevetr_id, albgrass_id, lai_id, \
        n_vegetation_region_urban, n_stream_sw_urban, n_stream_lw_urban, \
        sw_dn_direct_frac, air_ext_sw, air_ssa_sw, veg_ssa_sw, air_ext_lw, \
        air_ssa_lw, veg_ssa_lw, veg_fsd_const, veg_contact_fraction_const, \
        ground_albedo_dir_mult_fact, use_sw_direct_albedo, height, building_frac, \
        veg_frac, building_scale, veg_scale, alb_roof, emis_roof, alb_wall, \
        emis_wall, roof_albedo_dir_mult_fact, wall_specular_frac, alb_next, qn_surf, \
        qn_roof, qn_wall, qn_ind_snow, kup_ind_snow, tsurf_ind_snow, tsurf_ind, \
        dataoutlinespartacus)
    
    
    Defined at suews_ctrl_driver.fpp lines 1674-1866
    
    Parameters
    ----------
    storageheatmethod : int
    netradiationmethod : int
    snowuse : int
    tstep : int
    nlayer : int
    snowpack_prev : float array
    tau_a : float
    tau_f : float
    snowalbmax : float
    snowalbmin : float
    diagnose : int
    ldown_obs : float
    fcld_obs : float
    dectime : float
    zenith_deg : float
    tsurf_0 : float
    kdown : float
    tair_c : float
    avrh : float
    ea_hpa : float
    qn1_obs : float
    snowalb_prev : float
    snowfrac_prev : float array
    diagqn : int
    narp_trans_site : float
    narp_emis_snow : float
    icefrac : float array
    sfr_surf : float array
    sfr_roof : float array
    sfr_wall : float array
    tsfc_surf : float array
    tsfc_roof : float array
    tsfc_wall : float array
    emis : float array
    alb_prev : float array
    albdectr_id : float
    albevetr_id : float
    albgrass_id : float
    lai_id : float array
    n_vegetation_region_urban : int
    n_stream_sw_urban : int
    n_stream_lw_urban : int
    sw_dn_direct_frac : float
    air_ext_sw : float
    air_ssa_sw : float
    veg_ssa_sw : float
    air_ext_lw : float
    air_ssa_lw : float
    veg_ssa_lw : float
    veg_fsd_const : float
    veg_contact_fraction_const : float
    ground_albedo_dir_mult_fact : float
    use_sw_direct_albedo : bool
    height : float array
    building_frac : float array
    veg_frac : float array
    building_scale : float array
    veg_scale : float array
    alb_roof : float array
    emis_roof : float array
    alb_wall : float array
    emis_wall : float array
    roof_albedo_dir_mult_fact : float array
    wall_specular_frac : float array
    alb_next : float array
    qn_surf : float array
    qn_roof : float array
    qn_wall : float array
    qn_ind_snow : float array
    kup_ind_snow : float array
    tsurf_ind_snow : float array
    tsurf_ind : float array
    dataoutlinespartacus : float array
    
    Returns
    -------
    ldown : float
    fcld : float
    qn : float
    qn_snowfree : float
    qn_snow : float
    kclear : float
    kup : float
    lup : float
    tsurf : float
    albedo_snow : float
    snowalb_next : float
    
    """
    ldown, fcld, qn, qn_snowfree, qn_snow, kclear, kup, lup, tsurf, albedo_snow, \
        snowalb_next = \
        _suews_driver.f90wrap_suews_cal_qn(storageheatmethod=storageheatmethod, \
        netradiationmethod=netradiationmethod, snowuse=snowuse, tstep=tstep, \
        nlayer=nlayer, snowpack_prev=snowpack_prev, tau_a=tau_a, tau_f=tau_f, \
        snowalbmax=snowalbmax, snowalbmin=snowalbmin, diagnose=diagnose, \
        ldown_obs=ldown_obs, fcld_obs=fcld_obs, dectime=dectime, \
        zenith_deg=zenith_deg, tsurf_0=tsurf_0, kdown=kdown, tair_c=tair_c, \
        avrh=avrh, ea_hpa=ea_hpa, qn1_obs=qn1_obs, snowalb_prev=snowalb_prev, \
        snowfrac_prev=snowfrac_prev, diagqn=diagqn, narp_trans_site=narp_trans_site, \
        narp_emis_snow=narp_emis_snow, icefrac=icefrac, sfr_surf=sfr_surf, \
        sfr_roof=sfr_roof, sfr_wall=sfr_wall, tsfc_surf=tsfc_surf, \
        tsfc_roof=tsfc_roof, tsfc_wall=tsfc_wall, emis=emis, alb_prev=alb_prev, \
        albdectr_id=albdectr_id, albevetr_id=albevetr_id, albgrass_id=albgrass_id, \
        lai_id=lai_id, n_vegetation_region_urban=n_vegetation_region_urban, \
        n_stream_sw_urban=n_stream_sw_urban, n_stream_lw_urban=n_stream_lw_urban, \
        sw_dn_direct_frac=sw_dn_direct_frac, air_ext_sw=air_ext_sw, \
        air_ssa_sw=air_ssa_sw, veg_ssa_sw=veg_ssa_sw, air_ext_lw=air_ext_lw, \
        air_ssa_lw=air_ssa_lw, veg_ssa_lw=veg_ssa_lw, veg_fsd_const=veg_fsd_const, \
        veg_contact_fraction_const=veg_contact_fraction_const, \
        ground_albedo_dir_mult_fact=ground_albedo_dir_mult_fact, \
        use_sw_direct_albedo=use_sw_direct_albedo, height=height, \
        building_frac=building_frac, veg_frac=veg_frac, \
        building_scale=building_scale, veg_scale=veg_scale, alb_roof=alb_roof, \
        emis_roof=emis_roof, alb_wall=alb_wall, emis_wall=emis_wall, \
        roof_albedo_dir_mult_fact=roof_albedo_dir_mult_fact, \
        wall_specular_frac=wall_specular_frac, alb_next=alb_next, qn_surf=qn_surf, \
        qn_roof=qn_roof, qn_wall=qn_wall, qn_ind_snow=qn_ind_snow, \
        kup_ind_snow=kup_ind_snow, tsurf_ind_snow=tsurf_ind_snow, \
        tsurf_ind=tsurf_ind, dataoutlinespartacus=dataoutlinespartacus)
    return ldown, fcld, qn, qn_snowfree, qn_snow, kclear, kup, lup, tsurf, \
        albedo_snow, snowalb_next

def suews_cal_qs(storageheatmethod, qs_obs, ohmincqf, gridiv, id, tstep, \
    dt_since_start, diagnose, nlayer, qg_surf, qg_roof, qg_wall, tsfc_roof, \
    tin_roof, temp_in_roof, k_roof, cp_roof, dz_roof, sfr_roof, tsfc_wall, \
    tin_wall, temp_in_wall, k_wall, cp_wall, dz_wall, sfr_wall, tsfc_surf, \
    tin_surf, temp_in_surf, k_surf, cp_surf, dz_surf, sfr_surf, ohm_coef, \
    ohm_threshsw, ohm_threshwd, soilstore_id, soilstorecap, state_id, snowuse, \
    snowfrac, diagqs, hdd_id, metforcingdata_grid, ts5mindata_ir, qf, qn, avkdn, \
    avu1, temp_c, zenith_deg, avrh, press_hpa, ldown, bldgh, alb, emis, cpanohm, \
    kkanohm, chanohm, emissionsmethod, tair_av, qn_av_prev, dqndt_prev, \
    qn_s_av_prev, dqnsdt_prev, storedrainprm, dataoutlineestm, deltaqi, \
    temp_out_roof, qs_roof, temp_out_wall, qs_wall, temp_out_surf, qs_surf):
    """
    qn_s, qs, qn_av_next, dqndt_next, qn_s_av_next, dqnsdt_next, a1, a2, a3 = \
        suews_cal_qs(storageheatmethod, qs_obs, ohmincqf, gridiv, id, tstep, \
        dt_since_start, diagnose, nlayer, qg_surf, qg_roof, qg_wall, tsfc_roof, \
        tin_roof, temp_in_roof, k_roof, cp_roof, dz_roof, sfr_roof, tsfc_wall, \
        tin_wall, temp_in_wall, k_wall, cp_wall, dz_wall, sfr_wall, tsfc_surf, \
        tin_surf, temp_in_surf, k_surf, cp_surf, dz_surf, sfr_surf, ohm_coef, \
        ohm_threshsw, ohm_threshwd, soilstore_id, soilstorecap, state_id, snowuse, \
        snowfrac, diagqs, hdd_id, metforcingdata_grid, ts5mindata_ir, qf, qn, avkdn, \
        avu1, temp_c, zenith_deg, avrh, press_hpa, ldown, bldgh, alb, emis, cpanohm, \
        kkanohm, chanohm, emissionsmethod, tair_av, qn_av_prev, dqndt_prev, \
        qn_s_av_prev, dqnsdt_prev, storedrainprm, dataoutlineestm, deltaqi, \
        temp_out_roof, qs_roof, temp_out_wall, qs_wall, temp_out_surf, qs_surf)
    
    
    Defined at suews_ctrl_driver.fpp lines 1890-2084
    
    Parameters
    ----------
    storageheatmethod : int
    qs_obs : float
    ohmincqf : int
    gridiv : int
    id : int
    tstep : int
    dt_since_start : int
    diagnose : int
    nlayer : int
    qg_surf : float array
    qg_roof : float array
    qg_wall : float array
    tsfc_roof : float array
    tin_roof : float array
    temp_in_roof : float array
    k_roof : float array
    cp_roof : float array
    dz_roof : float array
    sfr_roof : float array
    tsfc_wall : float array
    tin_wall : float array
    temp_in_wall : float array
    k_wall : float array
    cp_wall : float array
    dz_wall : float array
    sfr_wall : float array
    tsfc_surf : float array
    tin_surf : float array
    temp_in_surf : float array
    k_surf : float array
    cp_surf : float array
    dz_surf : float array
    sfr_surf : float array
    ohm_coef : float array
    ohm_threshsw : float array
    ohm_threshwd : float array
    soilstore_id : float array
    soilstorecap : float array
    state_id : float array
    snowuse : int
    snowfrac : float array
    diagqs : int
    hdd_id : float array
    metforcingdata_grid : float array
    ts5mindata_ir : float array
    qf : float
    qn : float
    avkdn : float
    avu1 : float
    temp_c : float
    zenith_deg : float
    avrh : float
    press_hpa : float
    ldown : float
    bldgh : float
    alb : float array
    emis : float array
    cpanohm : float array
    kkanohm : float array
    chanohm : float array
    emissionsmethod : int
    tair_av : float
    qn_av_prev : float
    dqndt_prev : float
    qn_s_av_prev : float
    dqnsdt_prev : float
    storedrainprm : float array
    dataoutlineestm : float array
    deltaqi : float array
    temp_out_roof : float array
    qs_roof : float array
    temp_out_wall : float array
    qs_wall : float array
    temp_out_surf : float array
    qs_surf : float array
    
    Returns
    -------
    qn_s : float
    qs : float
    qn_av_next : float
    dqndt_next : float
    qn_s_av_next : float
    dqnsdt_next : float
    a1 : float
    a2 : float
    a3 : float
    
    """
    qn_s, qs, qn_av_next, dqndt_next, qn_s_av_next, dqnsdt_next, a1, a2, a3 = \
        _suews_driver.f90wrap_suews_cal_qs(storageheatmethod=storageheatmethod, \
        qs_obs=qs_obs, ohmincqf=ohmincqf, gridiv=gridiv, id=id, tstep=tstep, \
        dt_since_start=dt_since_start, diagnose=diagnose, nlayer=nlayer, \
        qg_surf=qg_surf, qg_roof=qg_roof, qg_wall=qg_wall, tsfc_roof=tsfc_roof, \
        tin_roof=tin_roof, temp_in_roof=temp_in_roof, k_roof=k_roof, \
        cp_roof=cp_roof, dz_roof=dz_roof, sfr_roof=sfr_roof, tsfc_wall=tsfc_wall, \
        tin_wall=tin_wall, temp_in_wall=temp_in_wall, k_wall=k_wall, \
        cp_wall=cp_wall, dz_wall=dz_wall, sfr_wall=sfr_wall, tsfc_surf=tsfc_surf, \
        tin_surf=tin_surf, temp_in_surf=temp_in_surf, k_surf=k_surf, \
        cp_surf=cp_surf, dz_surf=dz_surf, sfr_surf=sfr_surf, ohm_coef=ohm_coef, \
        ohm_threshsw=ohm_threshsw, ohm_threshwd=ohm_threshwd, \
        soilstore_id=soilstore_id, soilstorecap=soilstorecap, state_id=state_id, \
        snowuse=snowuse, snowfrac=snowfrac, diagqs=diagqs, hdd_id=hdd_id, \
        metforcingdata_grid=metforcingdata_grid, ts5mindata_ir=ts5mindata_ir, qf=qf, \
        qn=qn, avkdn=avkdn, avu1=avu1, temp_c=temp_c, zenith_deg=zenith_deg, \
        avrh=avrh, press_hpa=press_hpa, ldown=ldown, bldgh=bldgh, alb=alb, \
        emis=emis, cpanohm=cpanohm, kkanohm=kkanohm, chanohm=chanohm, \
        emissionsmethod=emissionsmethod, tair_av=tair_av, qn_av_prev=qn_av_prev, \
        dqndt_prev=dqndt_prev, qn_s_av_prev=qn_s_av_prev, dqnsdt_prev=dqnsdt_prev, \
        storedrainprm=storedrainprm, dataoutlineestm=dataoutlineestm, \
        deltaqi=deltaqi, temp_out_roof=temp_out_roof, qs_roof=qs_roof, \
        temp_out_wall=temp_out_wall, qs_wall=qs_wall, temp_out_surf=temp_out_surf, \
        qs_surf=qs_surf)
    return qn_s, qs, qn_av_next, dqndt_next, qn_s_av_next, dqnsdt_next, a1, a2, a3

def suews_cal_water(diagnose, snowuse, nonwaterfraction, addpipes, \
    addimpervious, addveg, addwaterbody, state_id, sfr_surf, storedrainprm, \
    waterdist, nsh_real, drain, frac_water2runoff, addwater):
    """
    drain_per_tstep, additionalwater, runoffpipes, runoff_per_interval = \
        suews_cal_water(diagnose, snowuse, nonwaterfraction, addpipes, \
        addimpervious, addveg, addwaterbody, state_id, sfr_surf, storedrainprm, \
        waterdist, nsh_real, drain, frac_water2runoff, addwater)
    
    
    Defined at suews_ctrl_driver.fpp lines 2095-2164
    
    Parameters
    ----------
    diagnose : int
    snowuse : int
    nonwaterfraction : float
    addpipes : float
    addimpervious : float
    addveg : float
    addwaterbody : float
    state_id : float array
    sfr_surf : float array
    storedrainprm : float array
    waterdist : float array
    nsh_real : float
    drain : float array
    frac_water2runoff : float array
    addwater : float array
    
    Returns
    -------
    drain_per_tstep : float
    additionalwater : float
    runoffpipes : float
    runoff_per_interval : float
    
    ============= Grid-to-grid runoff =============
     Calculate additional water coming from other grids
     i.e. the variables addImpervious, addVeg, addWaterBody, addPipes
    call RunoffFromGrid(GridFromFrac)
    Need to code between-grid water transfer
     Sum water coming from other grids(these are expressed as depths over the whole \
         surface)
    """
    drain_per_tstep, additionalwater, runoffpipes, runoff_per_interval = \
        _suews_driver.f90wrap_suews_cal_water(diagnose=diagnose, snowuse=snowuse, \
        nonwaterfraction=nonwaterfraction, addpipes=addpipes, \
        addimpervious=addimpervious, addveg=addveg, addwaterbody=addwaterbody, \
        state_id=state_id, sfr_surf=sfr_surf, storedrainprm=storedrainprm, \
        waterdist=waterdist, nsh_real=nsh_real, drain=drain, \
        frac_water2runoff=frac_water2runoff, addwater=addwater)
    return drain_per_tstep, additionalwater, runoffpipes, runoff_per_interval

def suews_init_qh(avdens, avcp, h_mod, qn1, dectime):
    """
    h_init = suews_init_qh(avdens, avcp, h_mod, qn1, dectime)
    
    
    Defined at suews_ctrl_driver.fpp lines 2170-2191
    
    Parameters
    ----------
    avdens : float
    avcp : float
    h_mod : float
    qn1 : float
    dectime : float
    
    Returns
    -------
    h_init : float
    
    """
    h_init = _suews_driver.f90wrap_suews_init_qh(avdens=avdens, avcp=avcp, \
        h_mod=h_mod, qn1=qn1, dectime=dectime)
    return h_init

def suews_cal_snow(diagnose, nlayer, tstep, imin, it, evapmethod, dayofweek_id, \
    crwmin, crwmax, dectime, avdens, avcp, lv_j_kg, lvs_j_kg, avrh, press_hpa, \
    temp_c, rasnow, psyc_hpa, sice_hpa, tau_r, radmeltfact, tempmeltfact, \
    snowalbmax, preciplimit, preciplimitalb, qn_ind_snow, kup_ind_snow, deltaqi, \
    tsurf_ind_snow, snowalb_in, pervfraction, vegfraction, addimpervious, \
    qn_snowfree, qf, qs, vpd_hpa, s_hpa, rs, ra, rb, snowdensmax, snowdensmin, \
    precip, pipecapacity, runofftowater, addveg, snowlimpaved, snowlimbldg, \
    flowchange, drain, wetthresh_surf, soilstorecap, tsurf_ind, sfr_surf, \
    addwater, addwaterrunoff, storedrainprm, snowpacklimit, snowprof_24hr, \
    snowpack_in, snowfrac_in, snowwater_in, icefrac_in, snowdens_in, \
    snowfallcum_in, state_id_in, soilstore_id_in, qn_surf, qs_surf, snowremoval, \
    snowpack_out, snowfrac_out, snowwater_out, icefrac_out, snowdens_out, \
    state_id_out, soilstore_id_out, qe_surf, qe_roof, qe_wall, rss_surf, \
    dataoutlinesnow):
    """
    snowfallcum_out, state_per_tstep, nwstate_per_tstep, qe, snowalb_out, swe, \
        chsnow_per_tstep, ev_per_tstep, runoff_per_tstep, surf_chang_per_tstep, \
        runoffpipes, mwstore, runoffwaterbody, runoffagveg, runoffagimpervious = \
        suews_cal_snow(diagnose, nlayer, tstep, imin, it, evapmethod, dayofweek_id, \
        crwmin, crwmax, dectime, avdens, avcp, lv_j_kg, lvs_j_kg, avrh, press_hpa, \
        temp_c, rasnow, psyc_hpa, sice_hpa, tau_r, radmeltfact, tempmeltfact, \
        snowalbmax, preciplimit, preciplimitalb, qn_ind_snow, kup_ind_snow, deltaqi, \
        tsurf_ind_snow, snowalb_in, pervfraction, vegfraction, addimpervious, \
        qn_snowfree, qf, qs, vpd_hpa, s_hpa, rs, ra, rb, snowdensmax, snowdensmin, \
        precip, pipecapacity, runofftowater, addveg, snowlimpaved, snowlimbldg, \
        flowchange, drain, wetthresh_surf, soilstorecap, tsurf_ind, sfr_surf, \
        addwater, addwaterrunoff, storedrainprm, snowpacklimit, snowprof_24hr, \
        snowpack_in, snowfrac_in, snowwater_in, icefrac_in, snowdens_in, \
        snowfallcum_in, state_id_in, soilstore_id_in, qn_surf, qs_surf, snowremoval, \
        snowpack_out, snowfrac_out, snowwater_out, icefrac_out, snowdens_out, \
        state_id_out, soilstore_id_out, qe_surf, qe_roof, qe_wall, rss_surf, \
        dataoutlinesnow)
    
    
    Defined at suews_ctrl_driver.fpp lines 2220-2521
    
    Parameters
    ----------
    diagnose : int
    nlayer : int
    tstep : int
    imin : int
    it : int
    evapmethod : int
    dayofweek_id : int array
    crwmin : float
    crwmax : float
    dectime : float
    avdens : float
    avcp : float
    lv_j_kg : float
    lvs_j_kg : float
    avrh : float
    press_hpa : float
    temp_c : float
    rasnow : float
    psyc_hpa : float
    sice_hpa : float
    tau_r : float
    radmeltfact : float
    tempmeltfact : float
    snowalbmax : float
    preciplimit : float
    preciplimitalb : float
    qn_ind_snow : float array
    kup_ind_snow : float array
    deltaqi : float array
    tsurf_ind_snow : float array
    snowalb_in : float
    pervfraction : float
    vegfraction : float
    addimpervious : float
    qn_snowfree : float
    qf : float
    qs : float
    vpd_hpa : float
    s_hpa : float
    rs : float
    ra : float
    rb : float
    snowdensmax : float
    snowdensmin : float
    precip : float
    pipecapacity : float
    runofftowater : float
    addveg : float
    snowlimpaved : float
    snowlimbldg : float
    flowchange : float
    drain : float array
    wetthresh_surf : float array
    soilstorecap : float array
    tsurf_ind : float array
    sfr_surf : float array
    addwater : float array
    addwaterrunoff : float array
    storedrainprm : float array
    snowpacklimit : float array
    snowprof_24hr : float array
    snowpack_in : float array
    snowfrac_in : float array
    snowwater_in : float array
    icefrac_in : float array
    snowdens_in : float array
    snowfallcum_in : float
    state_id_in : float array
    soilstore_id_in : float array
    qn_surf : float array
    qs_surf : float array
    snowremoval : float array
    snowpack_out : float array
    snowfrac_out : float array
    snowwater_out : float array
    icefrac_out : float array
    snowdens_out : float array
    state_id_out : float array
    soilstore_id_out : float array
    qe_surf : float array
    qe_roof : float array
    qe_wall : float array
    rss_surf : float array
    dataoutlinesnow : float array
    
    Returns
    -------
    snowfallcum_out : float
    state_per_tstep : float
    nwstate_per_tstep : float
    qe : float
    snowalb_out : float
    swe : float
    chsnow_per_tstep : float
    ev_per_tstep : float
    runoff_per_tstep : float
    surf_chang_per_tstep : float
    runoffpipes : float
    mwstore : float
    runoffwaterbody : float
    runoffagveg : float
    runoffagimpervious : float
    
    """
    snowfallcum_out, state_per_tstep, nwstate_per_tstep, qe, snowalb_out, swe, \
        chsnow_per_tstep, ev_per_tstep, runoff_per_tstep, surf_chang_per_tstep, \
        runoffpipes, mwstore, runoffwaterbody, runoffagveg, runoffagimpervious = \
        _suews_driver.f90wrap_suews_cal_snow(diagnose=diagnose, nlayer=nlayer, \
        tstep=tstep, imin=imin, it=it, evapmethod=evapmethod, \
        dayofweek_id=dayofweek_id, crwmin=crwmin, crwmax=crwmax, dectime=dectime, \
        avdens=avdens, avcp=avcp, lv_j_kg=lv_j_kg, lvs_j_kg=lvs_j_kg, avrh=avrh, \
        press_hpa=press_hpa, temp_c=temp_c, rasnow=rasnow, psyc_hpa=psyc_hpa, \
        sice_hpa=sice_hpa, tau_r=tau_r, radmeltfact=radmeltfact, \
        tempmeltfact=tempmeltfact, snowalbmax=snowalbmax, preciplimit=preciplimit, \
        preciplimitalb=preciplimitalb, qn_ind_snow=qn_ind_snow, \
        kup_ind_snow=kup_ind_snow, deltaqi=deltaqi, tsurf_ind_snow=tsurf_ind_snow, \
        snowalb_in=snowalb_in, pervfraction=pervfraction, vegfraction=vegfraction, \
        addimpervious=addimpervious, qn_snowfree=qn_snowfree, qf=qf, qs=qs, \
        vpd_hpa=vpd_hpa, s_hpa=s_hpa, rs=rs, ra=ra, rb=rb, snowdensmax=snowdensmax, \
        snowdensmin=snowdensmin, precip=precip, pipecapacity=pipecapacity, \
        runofftowater=runofftowater, addveg=addveg, snowlimpaved=snowlimpaved, \
        snowlimbldg=snowlimbldg, flowchange=flowchange, drain=drain, \
        wetthresh_surf=wetthresh_surf, soilstorecap=soilstorecap, \
        tsurf_ind=tsurf_ind, sfr_surf=sfr_surf, addwater=addwater, \
        addwaterrunoff=addwaterrunoff, storedrainprm=storedrainprm, \
        snowpacklimit=snowpacklimit, snowprof_24hr=snowprof_24hr, \
        snowpack_in=snowpack_in, snowfrac_in=snowfrac_in, snowwater_in=snowwater_in, \
        icefrac_in=icefrac_in, snowdens_in=snowdens_in, \
        snowfallcum_in=snowfallcum_in, state_id_in=state_id_in, \
        soilstore_id_in=soilstore_id_in, qn_surf=qn_surf, qs_surf=qs_surf, \
        snowremoval=snowremoval, snowpack_out=snowpack_out, \
        snowfrac_out=snowfrac_out, snowwater_out=snowwater_out, \
        icefrac_out=icefrac_out, snowdens_out=snowdens_out, \
        state_id_out=state_id_out, soilstore_id_out=soilstore_id_out, \
        qe_surf=qe_surf, qe_roof=qe_roof, qe_wall=qe_wall, rss_surf=rss_surf, \
        dataoutlinesnow=dataoutlinesnow)
    return snowfallcum_out, state_per_tstep, nwstate_per_tstep, qe, snowalb_out, \
        swe, chsnow_per_tstep, ev_per_tstep, runoff_per_tstep, surf_chang_per_tstep, \
        runoffpipes, mwstore, runoffwaterbody, runoffagveg, runoffagimpervious

def suews_cal_qe(diagnose, storageheatmethod, nlayer, tstep, evapmethod, avdens, \
    avcp, lv_j_kg, psyc_hpa, pervfraction, addimpervious, qf, vpd_hpa, s_hpa, \
    rs, ra_h, rb, precip, pipecapacity, runofftowater, nonwaterfraction, \
    wu_surf, addveg, addwaterbody, addwater_surf, flowchange, drain_surf, \
    frac_water2runoff_surf, storedrainprm, sfr_surf, statelimit_surf, \
    soilstorecap_surf, wetthresh_surf, state_surf_in, soilstore_surf_in, \
    qn_surf, qs_surf, sfr_roof, statelimit_roof, soilstorecap_roof, \
    wetthresh_roof, state_roof_in, soilstore_roof_in, qn_roof, qs_roof, \
    sfr_wall, statelimit_wall, soilstorecap_wall, wetthresh_wall, state_wall_in, \
    soilstore_wall_in, qn_wall, qs_wall, state_surf_out, soilstore_surf_out, \
    ev_surf, state_roof_out, soilstore_roof_out, ev_roof, state_wall_out, \
    soilstore_wall_out, ev_wall, ev0_surf, qe0_surf, qe_surf, qe_roof, qe_wall, \
    rss_surf):
    """
    state_grid, nwstate_grid, qe, ev_grid, runoff_grid, surf_chang_grid, \
        runoffpipes_grid, runoffwaterbody_grid, runoffagveg_grid, \
        runoffagimpervious_grid = suews_cal_qe(diagnose, storageheatmethod, nlayer, \
        tstep, evapmethod, avdens, avcp, lv_j_kg, psyc_hpa, pervfraction, \
        addimpervious, qf, vpd_hpa, s_hpa, rs, ra_h, rb, precip, pipecapacity, \
        runofftowater, nonwaterfraction, wu_surf, addveg, addwaterbody, \
        addwater_surf, flowchange, drain_surf, frac_water2runoff_surf, \
        storedrainprm, sfr_surf, statelimit_surf, soilstorecap_surf, wetthresh_surf, \
        state_surf_in, soilstore_surf_in, qn_surf, qs_surf, sfr_roof, \
        statelimit_roof, soilstorecap_roof, wetthresh_roof, state_roof_in, \
        soilstore_roof_in, qn_roof, qs_roof, sfr_wall, statelimit_wall, \
        soilstorecap_wall, wetthresh_wall, state_wall_in, soilstore_wall_in, \
        qn_wall, qs_wall, state_surf_out, soilstore_surf_out, ev_surf, \
        state_roof_out, soilstore_roof_out, ev_roof, state_wall_out, \
        soilstore_wall_out, ev_wall, ev0_surf, qe0_surf, qe_surf, qe_roof, qe_wall, \
        rss_surf)
    
    
    Defined at suews_ctrl_driver.fpp lines 2553-2812
    
    Parameters
    ----------
    diagnose : int
    storageheatmethod : int
    nlayer : int
    tstep : int
    evapmethod : int
    avdens : float
    avcp : float
    lv_j_kg : float
    psyc_hpa : float
    pervfraction : float
    addimpervious : float
    qf : float
    vpd_hpa : float
    s_hpa : float
    rs : float
    ra_h : float
    rb : float
    precip : float
    pipecapacity : float
    runofftowater : float
    nonwaterfraction : float
    wu_surf : float array
    addveg : float
    addwaterbody : float
    addwater_surf : float array
    flowchange : float
    drain_surf : float array
    frac_water2runoff_surf : float array
    storedrainprm : float array
    sfr_surf : float array
    statelimit_surf : float array
    soilstorecap_surf : float array
    wetthresh_surf : float array
    state_surf_in : float array
    soilstore_surf_in : float array
    qn_surf : float array
    qs_surf : float array
    sfr_roof : float array
    statelimit_roof : float array
    soilstorecap_roof : float array
    wetthresh_roof : float array
    state_roof_in : float array
    soilstore_roof_in : float array
    qn_roof : float array
    qs_roof : float array
    sfr_wall : float array
    statelimit_wall : float array
    soilstorecap_wall : float array
    wetthresh_wall : float array
    state_wall_in : float array
    soilstore_wall_in : float array
    qn_wall : float array
    qs_wall : float array
    state_surf_out : float array
    soilstore_surf_out : float array
    ev_surf : float array
    state_roof_out : float array
    soilstore_roof_out : float array
    ev_roof : float array
    state_wall_out : float array
    soilstore_wall_out : float array
    ev_wall : float array
    ev0_surf : float array
    qe0_surf : float array
    qe_surf : float array
    qe_roof : float array
    qe_wall : float array
    rss_surf : float array
    
    Returns
    -------
    state_grid : float
    nwstate_grid : float
    qe : float
    ev_grid : float
    runoff_grid : float
    surf_chang_grid : float
    runoffpipes_grid : float
    runoffwaterbody_grid : float
    runoffagveg_grid : float
    runoffagimpervious_grid : float
    
    """
    state_grid, nwstate_grid, qe, ev_grid, runoff_grid, surf_chang_grid, \
        runoffpipes_grid, runoffwaterbody_grid, runoffagveg_grid, \
        runoffagimpervious_grid = \
        _suews_driver.f90wrap_suews_cal_qe(diagnose=diagnose, \
        storageheatmethod=storageheatmethod, nlayer=nlayer, tstep=tstep, \
        evapmethod=evapmethod, avdens=avdens, avcp=avcp, lv_j_kg=lv_j_kg, \
        psyc_hpa=psyc_hpa, pervfraction=pervfraction, addimpervious=addimpervious, \
        qf=qf, vpd_hpa=vpd_hpa, s_hpa=s_hpa, rs=rs, ra_h=ra_h, rb=rb, precip=precip, \
        pipecapacity=pipecapacity, runofftowater=runofftowater, \
        nonwaterfraction=nonwaterfraction, wu_surf=wu_surf, addveg=addveg, \
        addwaterbody=addwaterbody, addwater_surf=addwater_surf, \
        flowchange=flowchange, drain_surf=drain_surf, \
        frac_water2runoff_surf=frac_water2runoff_surf, storedrainprm=storedrainprm, \
        sfr_surf=sfr_surf, statelimit_surf=statelimit_surf, \
        soilstorecap_surf=soilstorecap_surf, wetthresh_surf=wetthresh_surf, \
        state_surf_in=state_surf_in, soilstore_surf_in=soilstore_surf_in, \
        qn_surf=qn_surf, qs_surf=qs_surf, sfr_roof=sfr_roof, \
        statelimit_roof=statelimit_roof, soilstorecap_roof=soilstorecap_roof, \
        wetthresh_roof=wetthresh_roof, state_roof_in=state_roof_in, \
        soilstore_roof_in=soilstore_roof_in, qn_roof=qn_roof, qs_roof=qs_roof, \
        sfr_wall=sfr_wall, statelimit_wall=statelimit_wall, \
        soilstorecap_wall=soilstorecap_wall, wetthresh_wall=wetthresh_wall, \
        state_wall_in=state_wall_in, soilstore_wall_in=soilstore_wall_in, \
        qn_wall=qn_wall, qs_wall=qs_wall, state_surf_out=state_surf_out, \
        soilstore_surf_out=soilstore_surf_out, ev_surf=ev_surf, \
        state_roof_out=state_roof_out, soilstore_roof_out=soilstore_roof_out, \
        ev_roof=ev_roof, state_wall_out=state_wall_out, \
        soilstore_wall_out=soilstore_wall_out, ev_wall=ev_wall, ev0_surf=ev0_surf, \
        qe0_surf=qe0_surf, qe_surf=qe_surf, qe_roof=qe_roof, qe_wall=qe_wall, \
        rss_surf=rss_surf)
    return state_grid, nwstate_grid, qe, ev_grid, runoff_grid, surf_chang_grid, \
        runoffpipes_grid, runoffwaterbody_grid, runoffagveg_grid, \
        runoffagimpervious_grid

def suews_cal_qh(qhmethod, nlayer, storageheatmethod, qn, qf, qmrain, qe, qs, \
    qmfreez, qm, avdens, avcp, sfr_surf, sfr_roof, sfr_wall, tsfc_surf, \
    tsfc_roof, tsfc_wall, temp_c, ra, qh_resist_surf, qh_resist_roof, \
    qh_resist_wall):
    """
    qh, qh_residual, qh_resist = suews_cal_qh(qhmethod, nlayer, storageheatmethod, \
        qn, qf, qmrain, qe, qs, qmfreez, qm, avdens, avcp, sfr_surf, sfr_roof, \
        sfr_wall, tsfc_surf, tsfc_roof, tsfc_wall, temp_c, ra, qh_resist_surf, \
        qh_resist_roof, qh_resist_wall)
    
    
    Defined at suews_ctrl_driver.fpp lines 2824-2890
    
    Parameters
    ----------
    qhmethod : int
    nlayer : int
    storageheatmethod : int
    qn : float
    qf : float
    qmrain : float
    qe : float
    qs : float
    qmfreez : float
    qm : float
    avdens : float
    avcp : float
    sfr_surf : float array
    sfr_roof : float array
    sfr_wall : float array
    tsfc_surf : float array
    tsfc_roof : float array
    tsfc_wall : float array
    temp_c : float
    ra : float
    qh_resist_surf : float array
    qh_resist_roof : float array
    qh_resist_wall : float array
    
    Returns
    -------
    qh : float
    qh_residual : float
    qh_resist : float
    
    """
    qh, qh_residual, qh_resist = \
        _suews_driver.f90wrap_suews_cal_qh(qhmethod=qhmethod, nlayer=nlayer, \
        storageheatmethod=storageheatmethod, qn=qn, qf=qf, qmrain=qmrain, qe=qe, \
        qs=qs, qmfreez=qmfreez, qm=qm, avdens=avdens, avcp=avcp, sfr_surf=sfr_surf, \
        sfr_roof=sfr_roof, sfr_wall=sfr_wall, tsfc_surf=tsfc_surf, \
        tsfc_roof=tsfc_roof, tsfc_wall=tsfc_wall, temp_c=temp_c, ra=ra, \
        qh_resist_surf=qh_resist_surf, qh_resist_roof=qh_resist_roof, \
        qh_resist_wall=qh_resist_wall)
    return qh, qh_residual, qh_resist

def suews_cal_resistance(stabilitymethod, diagnose, aerodynamicresistancemethod, \
    roughlenheatmethod, snowuse, id, it, gsmodel, smdmethod, avdens, avcp, \
    qh_init, zzd, z0m, zdm, avu1, temp_c, vegfraction, avkdn, kmax, g_max, g_k, \
    g_q_base, g_q_shape, g_t, g_sm, s1, s2, th, tl, dq, xsmd, vsmd, \
    maxconductance, laimax, lai_id, snowfrac, sfr_surf):
    """
    ustar, tstar, l_mod, zl, gsc, rs, ra, rasnow, rb, z0v, z0vsnow = \
        suews_cal_resistance(stabilitymethod, diagnose, aerodynamicresistancemethod, \
        roughlenheatmethod, snowuse, id, it, gsmodel, smdmethod, avdens, avcp, \
        qh_init, zzd, z0m, zdm, avu1, temp_c, vegfraction, avkdn, kmax, g_max, g_k, \
        g_q_base, g_q_shape, g_t, g_sm, s1, s2, th, tl, dq, xsmd, vsmd, \
        maxconductance, laimax, lai_id, snowfrac, sfr_surf)
    
    
    Defined at suews_ctrl_driver.fpp lines 2906-3021
    
    Parameters
    ----------
    stabilitymethod : int
    diagnose : int
    aerodynamicresistancemethod : int
    roughlenheatmethod : int
    snowuse : int
    id : int
    it : int
    gsmodel : int
    smdmethod : int
    avdens : float
    avcp : float
    qh_init : float
    zzd : float
    z0m : float
    zdm : float
    avu1 : float
    temp_c : float
    vegfraction : float
    avkdn : float
    kmax : float
    g_max : float
    g_k : float
    g_q_base : float
    g_q_shape : float
    g_t : float
    g_sm : float
    s1 : float
    s2 : float
    th : float
    tl : float
    dq : float
    xsmd : float
    vsmd : float
    maxconductance : float array
    laimax : float array
    lai_id : float array
    snowfrac : float array
    sfr_surf : float array
    
    Returns
    -------
    ustar : float
    tstar : float
    l_mod : float
    zl : float
    gsc : float
    rs : float
    ra : float
    rasnow : float
    rb : float
    z0v : float
    z0vsnow : float
    
    """
    ustar, tstar, l_mod, zl, gsc, rs, ra, rasnow, rb, z0v, z0vsnow = \
        _suews_driver.f90wrap_suews_cal_resistance(stabilitymethod=stabilitymethod, \
        diagnose=diagnose, aerodynamicresistancemethod=aerodynamicresistancemethod, \
        roughlenheatmethod=roughlenheatmethod, snowuse=snowuse, id=id, it=it, \
        gsmodel=gsmodel, smdmethod=smdmethod, avdens=avdens, avcp=avcp, \
        qh_init=qh_init, zzd=zzd, z0m=z0m, zdm=zdm, avu1=avu1, temp_c=temp_c, \
        vegfraction=vegfraction, avkdn=avkdn, kmax=kmax, g_max=g_max, g_k=g_k, \
        g_q_base=g_q_base, g_q_shape=g_q_shape, g_t=g_t, g_sm=g_sm, s1=s1, s2=s2, \
        th=th, tl=tl, dq=dq, xsmd=xsmd, vsmd=vsmd, maxconductance=maxconductance, \
        laimax=laimax, lai_id=lai_id, snowfrac=snowfrac, sfr_surf=sfr_surf)
    return ustar, tstar, l_mod, zl, gsc, rs, ra, rasnow, rb, z0v, z0vsnow

def suews_update_outputline(additionalwater, alb, avkdn, avu10_ms, azimuth, \
    chsnow_per_interval, dectime, drain_per_tstep, e_mod, ev_per_tstep, ext_wu, \
    fc, fc_build, fcld, fc_metab, fc_photo, fc_respi, fc_point, fc_traff, \
    flowchange, h_mod, id, imin, int_wu, it, iy, kup, lai_id, ldown, l_mod, lup, \
    mwh, mwstore, nsh_real, nwstate_per_tstep, precip, q2_gkg, qeout, qf, qh, \
    qh_resist, qm, qmfreez, qmrain, qn, qn_snow, qn_snowfree, qs, ra, \
    resistsurf, rh2, runoffagimpervious, runoffagveg, runoff_per_tstep, \
    runoffpipes, runoffsoil_per_tstep, runoffwaterbody, sfr_surf, smd, \
    smd_nsurf, snowalb, snowremoval, state_id, state_per_tstep, \
    surf_chang_per_tstep, swe, t2_c, tskin_c, tot_chang_per_tstep, tsurf, ustar, \
    wu_nsurf, z0m, zdm, zenith_deg, datetimeline, dataoutlinesuews):
    """
    suews_update_outputline(additionalwater, alb, avkdn, avu10_ms, azimuth, \
        chsnow_per_interval, dectime, drain_per_tstep, e_mod, ev_per_tstep, ext_wu, \
        fc, fc_build, fcld, fc_metab, fc_photo, fc_respi, fc_point, fc_traff, \
        flowchange, h_mod, id, imin, int_wu, it, iy, kup, lai_id, ldown, l_mod, lup, \
        mwh, mwstore, nsh_real, nwstate_per_tstep, precip, q2_gkg, qeout, qf, qh, \
        qh_resist, qm, qmfreez, qmrain, qn, qn_snow, qn_snowfree, qs, ra, \
        resistsurf, rh2, runoffagimpervious, runoffagveg, runoff_per_tstep, \
        runoffpipes, runoffsoil_per_tstep, runoffwaterbody, sfr_surf, smd, \
        smd_nsurf, snowalb, snowremoval, state_id, state_per_tstep, \
        surf_chang_per_tstep, swe, t2_c, tskin_c, tot_chang_per_tstep, tsurf, ustar, \
        wu_nsurf, z0m, zdm, zenith_deg, datetimeline, dataoutlinesuews)
    
    
    Defined at suews_ctrl_driver.fpp lines 3043-3191
    
    Parameters
    ----------
    additionalwater : float
    alb : float array
    avkdn : float
    avu10_ms : float
    azimuth : float
    chsnow_per_interval : float
    dectime : float
    drain_per_tstep : float
    e_mod : float
    ev_per_tstep : float
    ext_wu : float
    fc : float
    fc_build : float
    fcld : float
    fc_metab : float
    fc_photo : float
    fc_respi : float
    fc_point : float
    fc_traff : float
    flowchange : float
    h_mod : float
    id : int
    imin : int
    int_wu : float
    it : int
    iy : int
    kup : float
    lai_id : float array
    ldown : float
    l_mod : float
    lup : float
    mwh : float
    mwstore : float
    nsh_real : float
    nwstate_per_tstep : float
    precip : float
    q2_gkg : float
    qeout : float
    qf : float
    qh : float
    qh_resist : float
    qm : float
    qmfreez : float
    qmrain : float
    qn : float
    qn_snow : float
    qn_snowfree : float
    qs : float
    ra : float
    resistsurf : float
    rh2 : float
    runoffagimpervious : float
    runoffagveg : float
    runoff_per_tstep : float
    runoffpipes : float
    runoffsoil_per_tstep : float
    runoffwaterbody : float
    sfr_surf : float array
    smd : float
    smd_nsurf : float array
    snowalb : float
    snowremoval : float array
    state_id : float array
    state_per_tstep : float
    surf_chang_per_tstep : float
    swe : float
    t2_c : float
    tskin_c : float
    tot_chang_per_tstep : float
    tsurf : float
    ustar : float
    wu_nsurf : float array
    z0m : float
    zdm : float
    zenith_deg : float
    datetimeline : float array
    dataoutlinesuews : float array
    
    =====================================================================
    ====================== Prepare data for output ======================
     values outside of reasonable range are set as NAN-like numbers. TS 10 Jun 2018
     Remove non-existing surface type from surface and soil outputs
     Added back in with NANs by HCW 24 Aug 2016
    """
    _suews_driver.f90wrap_suews_update_outputline(additionalwater=additionalwater, \
        alb=alb, avkdn=avkdn, avu10_ms=avu10_ms, azimuth=azimuth, \
        chsnow_per_interval=chsnow_per_interval, dectime=dectime, \
        drain_per_tstep=drain_per_tstep, e_mod=e_mod, ev_per_tstep=ev_per_tstep, \
        ext_wu=ext_wu, fc=fc, fc_build=fc_build, fcld=fcld, fc_metab=fc_metab, \
        fc_photo=fc_photo, fc_respi=fc_respi, fc_point=fc_point, fc_traff=fc_traff, \
        flowchange=flowchange, h_mod=h_mod, id=id, imin=imin, int_wu=int_wu, it=it, \
        iy=iy, kup=kup, lai_id=lai_id, ldown=ldown, l_mod=l_mod, lup=lup, mwh=mwh, \
        mwstore=mwstore, nsh_real=nsh_real, nwstate_per_tstep=nwstate_per_tstep, \
        precip=precip, q2_gkg=q2_gkg, qeout=qeout, qf=qf, qh=qh, \
        qh_resist=qh_resist, qm=qm, qmfreez=qmfreez, qmrain=qmrain, qn=qn, \
        qn_snow=qn_snow, qn_snowfree=qn_snowfree, qs=qs, ra=ra, \
        resistsurf=resistsurf, rh2=rh2, runoffagimpervious=runoffagimpervious, \
        runoffagveg=runoffagveg, runoff_per_tstep=runoff_per_tstep, \
        runoffpipes=runoffpipes, runoffsoil_per_tstep=runoffsoil_per_tstep, \
        runoffwaterbody=runoffwaterbody, sfr_surf=sfr_surf, smd=smd, \
        smd_nsurf=smd_nsurf, snowalb=snowalb, snowremoval=snowremoval, \
        state_id=state_id, state_per_tstep=state_per_tstep, \
        surf_chang_per_tstep=surf_chang_per_tstep, swe=swe, t2_c=t2_c, \
        tskin_c=tskin_c, tot_chang_per_tstep=tot_chang_per_tstep, tsurf=tsurf, \
        ustar=ustar, wu_nsurf=wu_nsurf, z0m=z0m, zdm=zdm, zenith_deg=zenith_deg, \
        datetimeline=datetimeline, dataoutlinesuews=dataoutlinesuews)

def estmext_update_outputline(iy, id, it, imin, dectime, nlayer, tsfc_out_surf, \
    qs_surf, tsfc_out_roof, qn_roof, qs_roof, qe_roof, qh_roof, state_roof, \
    soilstore_roof, tsfc_out_wall, qn_wall, qs_wall, qe_wall, qh_wall, \
    state_wall, soilstore_wall, datetimeline, dataoutlineestmext):
    """
    estmext_update_outputline(iy, id, it, imin, dectime, nlayer, tsfc_out_surf, \
        qs_surf, tsfc_out_roof, qn_roof, qs_roof, qe_roof, qh_roof, state_roof, \
        soilstore_roof, tsfc_out_wall, qn_wall, qs_wall, qe_wall, qh_wall, \
        state_wall, soilstore_wall, datetimeline, dataoutlineestmext)
    
    
    Defined at suews_ctrl_driver.fpp lines 3212-3279
    
    Parameters
    ----------
    iy : int
    id : int
    it : int
    imin : int
    dectime : float
    nlayer : int
    tsfc_out_surf : float array
    qs_surf : float array
    tsfc_out_roof : float array
    qn_roof : float array
    qs_roof : float array
    qe_roof : float array
    qh_roof : float array
    state_roof : float array
    soilstore_roof : float array
    tsfc_out_wall : float array
    qn_wall : float array
    qs_wall : float array
    qe_wall : float array
    qh_wall : float array
    state_wall : float array
    soilstore_wall : float array
    datetimeline : float array
    dataoutlineestmext : float array
    
    ====================update output line end==============================
    """
    _suews_driver.f90wrap_estmext_update_outputline(iy=iy, id=id, it=it, imin=imin, \
        dectime=dectime, nlayer=nlayer, tsfc_out_surf=tsfc_out_surf, \
        qs_surf=qs_surf, tsfc_out_roof=tsfc_out_roof, qn_roof=qn_roof, \
        qs_roof=qs_roof, qe_roof=qe_roof, qh_roof=qh_roof, state_roof=state_roof, \
        soilstore_roof=soilstore_roof, tsfc_out_wall=tsfc_out_wall, qn_wall=qn_wall, \
        qs_wall=qs_wall, qe_wall=qe_wall, qh_wall=qh_wall, state_wall=state_wall, \
        soilstore_wall=soilstore_wall, datetimeline=datetimeline, \
        dataoutlineestmext=dataoutlineestmext)

def fill_result_x(res_valid, n_fill):
    """
    res_filled = fill_result_x(res_valid, n_fill)
    
    
    Defined at suews_ctrl_driver.fpp lines 3282-3288
    
    Parameters
    ----------
    res_valid : float array
    n_fill : int
    
    Returns
    -------
    res_filled : float array
    
    """
    res_filled = _suews_driver.f90wrap_fill_result_x(res_valid=res_valid, \
        n_fill=n_fill)
    return res_filled

def suews_update_output(snowuse, storageheatmethod, readlinesmetdata, \
    numberofgrids, ir, gridiv, datetimeline, dataoutlinesuews, dataoutlinesnow, \
    dataoutlineestm, dataoutlinersl, dataoutlinebeers, dataoutlinedebug, \
    dataoutlinespartacus, dataoutlineestmext, dataoutsuews, dataoutsnow, \
    dataoutestm, dataoutrsl, dataoutbeers, dataoutdebug, dataoutspartacus, \
    dataoutestmext):
    """
    suews_update_output(snowuse, storageheatmethod, readlinesmetdata, numberofgrids, \
        ir, gridiv, datetimeline, dataoutlinesuews, dataoutlinesnow, \
        dataoutlineestm, dataoutlinersl, dataoutlinebeers, dataoutlinedebug, \
        dataoutlinespartacus, dataoutlineestmext, dataoutsuews, dataoutsnow, \
        dataoutestm, dataoutrsl, dataoutbeers, dataoutdebug, dataoutspartacus, \
        dataoutestmext)
    
    
    Defined at suews_ctrl_driver.fpp lines 3298-3343
    
    Parameters
    ----------
    snowuse : int
    storageheatmethod : int
    readlinesmetdata : int
    numberofgrids : int
    ir : int
    gridiv : int
    datetimeline : float array
    dataoutlinesuews : float array
    dataoutlinesnow : float array
    dataoutlineestm : float array
    dataoutlinersl : float array
    dataoutlinebeers : float array
    dataoutlinedebug : float array
    dataoutlinespartacus : float array
    dataoutlineestmext : float array
    dataoutsuews : float array
    dataoutsnow : float array
    dataoutestm : float array
    dataoutrsl : float array
    dataoutbeers : float array
    dataoutdebug : float array
    dataoutspartacus : float array
    dataoutestmext : float array
    
    ====================== update output arrays ==============================
    Define the overall output matrix to be printed out step by step
    """
    _suews_driver.f90wrap_suews_update_output(snowuse=snowuse, \
        storageheatmethod=storageheatmethod, readlinesmetdata=readlinesmetdata, \
        numberofgrids=numberofgrids, ir=ir, gridiv=gridiv, \
        datetimeline=datetimeline, dataoutlinesuews=dataoutlinesuews, \
        dataoutlinesnow=dataoutlinesnow, dataoutlineestm=dataoutlineestm, \
        dataoutlinersl=dataoutlinersl, dataoutlinebeers=dataoutlinebeers, \
        dataoutlinedebug=dataoutlinedebug, \
        dataoutlinespartacus=dataoutlinespartacus, \
        dataoutlineestmext=dataoutlineestmext, dataoutsuews=dataoutsuews, \
        dataoutsnow=dataoutsnow, dataoutestm=dataoutestm, dataoutrsl=dataoutrsl, \
        dataoutbeers=dataoutbeers, dataoutdebug=dataoutdebug, \
        dataoutspartacus=dataoutspartacus, dataoutestmext=dataoutestmext)

def suews_cal_surf(storageheatmethod, netradiationmethod, nlayer, sfr_surf, \
    building_frac, building_scale, height, sfr_roof, sfr_wall):
    """
    vegfraction, impervfraction, pervfraction, nonwaterfraction = \
        suews_cal_surf(storageheatmethod, netradiationmethod, nlayer, sfr_surf, \
        building_frac, building_scale, height, sfr_roof, sfr_wall)
    
    
    Defined at suews_ctrl_driver.fpp lines 3413-3456
    
    Parameters
    ----------
    storageheatmethod : int
    netradiationmethod : int
    nlayer : int
    sfr_surf : float array
    building_frac : float array
    building_scale : float array
    height : float array
    sfr_roof : float array
    sfr_wall : float array
    
    Returns
    -------
    vegfraction : float
    impervfraction : float
    pervfraction : float
    nonwaterfraction : float
    
    """
    vegfraction, impervfraction, pervfraction, nonwaterfraction = \
        _suews_driver.f90wrap_suews_cal_surf(storageheatmethod=storageheatmethod, \
        netradiationmethod=netradiationmethod, nlayer=nlayer, sfr_surf=sfr_surf, \
        building_frac=building_frac, building_scale=building_scale, height=height, \
        sfr_roof=sfr_roof, sfr_wall=sfr_wall)
    return vegfraction, impervfraction, pervfraction, nonwaterfraction

def set_nan(x):
    """
    xx = set_nan(x)
    
    
    Defined at suews_ctrl_driver.fpp lines 3575-3588
    
    Parameters
    ----------
    x : float
    
    Returns
    -------
    xx : float
    
    """
    xx = _suews_driver.f90wrap_set_nan(x=x)
    return xx

def square(x):
    """
    xx = square(x)
    
    
    Defined at suews_ctrl_driver.fpp lines 3592-3599
    
    Parameters
    ----------
    x : float
    
    Returns
    -------
    xx : float
    
    """
    xx = _suews_driver.f90wrap_square(x=x)
    return xx

def square_real(x):
    """
    xx = square_real(x)
    
    
    Defined at suews_ctrl_driver.fpp lines 3601-3608
    
    Parameters
    ----------
    x : float
    
    Returns
    -------
    xx : float
    
    """
    xx = _suews_driver.f90wrap_square_real(x=x)
    return xx

def output_name_n(i):
    """
    name, group, aggreg, outlevel = output_name_n(i)
    
    
    Defined at suews_ctrl_driver.fpp lines 3610-3630
    
    Parameters
    ----------
    i : int
    
    Returns
    -------
    name : str
    group : str
    aggreg : str
    outlevel : int
    
    """
    name, group, aggreg, outlevel = _suews_driver.f90wrap_output_name_n(i=i)
    return name, group, aggreg, outlevel

def output_size():
    """
    nvar = output_size()
    
    
    Defined at suews_ctrl_driver.fpp lines 3632-3703
    
    
    Returns
    -------
    nvar : int
    
    """
    nvar = _suews_driver.f90wrap_output_size()
    return nvar

def suews_cal_multitsteps(metforcingblock, len_sim, ah_min, ahprof_24hr, \
    ah_slope_cooling, ah_slope_heating, alb, albmax_dectr, albmax_evetr, \
    albmax_grass, albmin_dectr, albmin_evetr, albmin_grass, alpha_bioco2, \
    alpha_enh_bioco2, alt, baset, basete, beta_bioco2, beta_enh_bioco2, bldgh, \
    capmax_dec, capmin_dec, chanohm, co2pointsource, cpanohm, crwmax, crwmin, \
    daywat, daywatper, dectreeh, diagmethod, diagnose, drainrt, dt_since_start, \
    dqndt, qn_av, dqnsdt, qn_s_av, ef_umolco2perj, emis, emissionsmethod, \
    enef_v_jkm, enddls, evetreeh, faibldg, faidectree, faievetree, faut, \
    fcef_v_kgkm, flowchange, frfossilfuel_heat, frfossilfuel_nonheat, g_max, \
    g_k, g_q_base, g_q_shape, g_t, g_sm, gdd_id, gddfull, gridiv, gsmodel, \
    h_maintain, hdd_id, humactivity_24hr, icefrac, ie_a, ie_end, ie_m, ie_start, \
    internalwateruse_h, irrfracpaved, irrfracbldgs, irrfracevetr, irrfracdectr, \
    irrfracgrass, irrfracbsoil, irrfracwater, kkanohm, kmax, lai_id, laimax, \
    laimin, laipower, laitype, lat, lng, maxconductance, maxfcmetab, maxqfmetab, \
    snowwater, minfcmetab, minqfmetab, min_res_bioco2, narp_emis_snow, \
    narp_trans_site, netradiationmethod, ohm_coef, ohmincqf, ohm_threshsw, \
    ohm_threshwd, pipecapacity, popdensdaytime, popdensnighttime, popprof_24hr, \
    pormax_dec, pormin_dec, preciplimit, preciplimitalb, qf0_beu, qf_a, qf_b, \
    qf_c, nlayer, n_vegetation_region_urban, n_stream_sw_urban, \
    n_stream_lw_urban, sw_dn_direct_frac, air_ext_sw, air_ssa_sw, veg_ssa_sw, \
    air_ext_lw, air_ssa_lw, veg_ssa_lw, veg_fsd_const, \
    veg_contact_fraction_const, ground_albedo_dir_mult_fact, \
    use_sw_direct_albedo, height, building_frac, veg_frac, building_scale, \
    veg_scale, alb_roof, emis_roof, alb_wall, emis_wall, \
    roof_albedo_dir_mult_fact, wall_specular_frac, radmeltfact, raincover, \
    rainmaxres, resp_a, resp_b, roughlenheatmethod, roughlenmommethod, \
    runofftowater, s1, s2, sathydraulicconduct, sddfull, sdd_id, smdmethod, \
    snowalb, snowalbmax, snowalbmin, snowpacklimit, snowdens, snowdensmax, \
    snowdensmin, snowfallcum, snowfrac, snowlimbldg, snowlimpaved, snowpack, \
    snowprof_24hr, snowuse, soildepth, stabilitymethod, startdls, \
    soilstore_surf, soilstorecap_surf, state_surf, statelimit_surf, \
    wetthresh_surf, soilstore_roof, soilstorecap_roof, state_roof, \
    statelimit_roof, wetthresh_roof, soilstore_wall, soilstorecap_wall, \
    state_wall, statelimit_wall, wetthresh_wall, storageheatmethod, \
    storedrainprm, surfacearea, tair_av, tau_a, tau_f, tau_r, baset_cooling, \
    baset_heating, tempmeltfact, th, theta_bioco2, timezone, tl, trafficrate, \
    trafficunits, sfr_surf, tsfc_roof, tsfc_wall, tsfc_surf, temp_roof, \
    temp_wall, temp_surf, tin_roof, tin_wall, tin_surf, k_wall, k_roof, k_surf, \
    cp_wall, cp_roof, cp_surf, dz_wall, dz_roof, dz_surf, tmin_id, tmax_id, \
    lenday_id, traffprof_24hr, ts5mindata_ir, tstep, tstep_prev, veg_type, \
    waterdist, waterusemethod, wuday_id, decidcap_id, albdectr_id, albevetr_id, \
    albgrass_id, porosity_id, wuprofa_24hr, wuprofm_24hr, z, z0m_in, zdm_in, \
    dataoutblocksuews, dataoutblocksnow, dataoutblockestm, dataoutblockrsl, \
    dataoutblockbeers, dataoutblockdebug, dataoutblockspartacus, \
    dataoutblockestmext, dailystateblock):
    """
    suews_cal_multitsteps(metforcingblock, len_sim, ah_min, ahprof_24hr, \
        ah_slope_cooling, ah_slope_heating, alb, albmax_dectr, albmax_evetr, \
        albmax_grass, albmin_dectr, albmin_evetr, albmin_grass, alpha_bioco2, \
        alpha_enh_bioco2, alt, baset, basete, beta_bioco2, beta_enh_bioco2, bldgh, \
        capmax_dec, capmin_dec, chanohm, co2pointsource, cpanohm, crwmax, crwmin, \
        daywat, daywatper, dectreeh, diagmethod, diagnose, drainrt, dt_since_start, \
        dqndt, qn_av, dqnsdt, qn_s_av, ef_umolco2perj, emis, emissionsmethod, \
        enef_v_jkm, enddls, evetreeh, faibldg, faidectree, faievetree, faut, \
        fcef_v_kgkm, flowchange, frfossilfuel_heat, frfossilfuel_nonheat, g_max, \
        g_k, g_q_base, g_q_shape, g_t, g_sm, gdd_id, gddfull, gridiv, gsmodel, \
        h_maintain, hdd_id, humactivity_24hr, icefrac, ie_a, ie_end, ie_m, ie_start, \
        internalwateruse_h, irrfracpaved, irrfracbldgs, irrfracevetr, irrfracdectr, \
        irrfracgrass, irrfracbsoil, irrfracwater, kkanohm, kmax, lai_id, laimax, \
        laimin, laipower, laitype, lat, lng, maxconductance, maxfcmetab, maxqfmetab, \
        snowwater, minfcmetab, minqfmetab, min_res_bioco2, narp_emis_snow, \
        narp_trans_site, netradiationmethod, ohm_coef, ohmincqf, ohm_threshsw, \
        ohm_threshwd, pipecapacity, popdensdaytime, popdensnighttime, popprof_24hr, \
        pormax_dec, pormin_dec, preciplimit, preciplimitalb, qf0_beu, qf_a, qf_b, \
        qf_c, nlayer, n_vegetation_region_urban, n_stream_sw_urban, \
        n_stream_lw_urban, sw_dn_direct_frac, air_ext_sw, air_ssa_sw, veg_ssa_sw, \
        air_ext_lw, air_ssa_lw, veg_ssa_lw, veg_fsd_const, \
        veg_contact_fraction_const, ground_albedo_dir_mult_fact, \
        use_sw_direct_albedo, height, building_frac, veg_frac, building_scale, \
        veg_scale, alb_roof, emis_roof, alb_wall, emis_wall, \
        roof_albedo_dir_mult_fact, wall_specular_frac, radmeltfact, raincover, \
        rainmaxres, resp_a, resp_b, roughlenheatmethod, roughlenmommethod, \
        runofftowater, s1, s2, sathydraulicconduct, sddfull, sdd_id, smdmethod, \
        snowalb, snowalbmax, snowalbmin, snowpacklimit, snowdens, snowdensmax, \
        snowdensmin, snowfallcum, snowfrac, snowlimbldg, snowlimpaved, snowpack, \
        snowprof_24hr, snowuse, soildepth, stabilitymethod, startdls, \
        soilstore_surf, soilstorecap_surf, state_surf, statelimit_surf, \
        wetthresh_surf, soilstore_roof, soilstorecap_roof, state_roof, \
        statelimit_roof, wetthresh_roof, soilstore_wall, soilstorecap_wall, \
        state_wall, statelimit_wall, wetthresh_wall, storageheatmethod, \
        storedrainprm, surfacearea, tair_av, tau_a, tau_f, tau_r, baset_cooling, \
        baset_heating, tempmeltfact, th, theta_bioco2, timezone, tl, trafficrate, \
        trafficunits, sfr_surf, tsfc_roof, tsfc_wall, tsfc_surf, temp_roof, \
        temp_wall, temp_surf, tin_roof, tin_wall, tin_surf, k_wall, k_roof, k_surf, \
        cp_wall, cp_roof, cp_surf, dz_wall, dz_roof, dz_surf, tmin_id, tmax_id, \
        lenday_id, traffprof_24hr, ts5mindata_ir, tstep, tstep_prev, veg_type, \
        waterdist, waterusemethod, wuday_id, decidcap_id, albdectr_id, albevetr_id, \
        albgrass_id, porosity_id, wuprofa_24hr, wuprofm_24hr, z, z0m_in, zdm_in, \
        dataoutblocksuews, dataoutblocksnow, dataoutblockestm, dataoutblockrsl, \
        dataoutblockbeers, dataoutblockdebug, dataoutblockspartacus, \
        dataoutblockestmext, dailystateblock)
    
    
    Defined at suews_ctrl_driver.fpp lines 3705-4408
    
    Parameters
    ----------
    metforcingblock : float array
    len_sim : int
    ah_min : float array
    ahprof_24hr : float array
    ah_slope_cooling : float array
    ah_slope_heating : float array
    alb : float array
    albmax_dectr : float
    albmax_evetr : float
    albmax_grass : float
    albmin_dectr : float
    albmin_evetr : float
    albmin_grass : float
    alpha_bioco2 : float array
    alpha_enh_bioco2 : float array
    alt : float
    baset : float array
    basete : float array
    beta_bioco2 : float array
    beta_enh_bioco2 : float array
    bldgh : float
    capmax_dec : float
    capmin_dec : float
    chanohm : float array
    co2pointsource : float
    cpanohm : float array
    crwmax : float
    crwmin : float
    daywat : float array
    daywatper : float array
    dectreeh : float
    diagmethod : int
    diagnose : int
    drainrt : float
    dt_since_start : int
    dqndt : float
    qn_av : float
    dqnsdt : float
    qn_s_av : float
    ef_umolco2perj : float
    emis : float array
    emissionsmethod : int
    enef_v_jkm : float
    enddls : int
    evetreeh : float
    faibldg : float
    faidectree : float
    faievetree : float
    faut : float
    fcef_v_kgkm : float array
    flowchange : float
    frfossilfuel_heat : float
    frfossilfuel_nonheat : float
    g_max : float
    g_k : float
    g_q_base : float
    g_q_shape : float
    g_t : float
    g_sm : float
    gdd_id : float array
    gddfull : float array
    gridiv : int
    gsmodel : int
    h_maintain : float
    hdd_id : float array
    humactivity_24hr : float array
    icefrac : float array
    ie_a : float array
    ie_end : int
    ie_m : float array
    ie_start : int
    internalwateruse_h : float
    irrfracpaved : float
    irrfracbldgs : float
    irrfracevetr : float
    irrfracdectr : float
    irrfracgrass : float
    irrfracbsoil : float
    irrfracwater : float
    kkanohm : float array
    kmax : float
    lai_id : float array
    laimax : float array
    laimin : float array
    laipower : float array
    laitype : int array
    lat : float
    lng : float
    maxconductance : float array
    maxfcmetab : float
    maxqfmetab : float
    snowwater : float array
    minfcmetab : float
    minqfmetab : float
    min_res_bioco2 : float array
    narp_emis_snow : float
    narp_trans_site : float
    netradiationmethod : int
    ohm_coef : float array
    ohmincqf : int
    ohm_threshsw : float array
    ohm_threshwd : float array
    pipecapacity : float
    popdensdaytime : float array
    popdensnighttime : float
    popprof_24hr : float array
    pormax_dec : float
    pormin_dec : float
    preciplimit : float
    preciplimitalb : float
    qf0_beu : float array
    qf_a : float array
    qf_b : float array
    qf_c : float array
    nlayer : int
    n_vegetation_region_urban : int
    n_stream_sw_urban : int
    n_stream_lw_urban : int
    sw_dn_direct_frac : float
    air_ext_sw : float
    air_ssa_sw : float
    veg_ssa_sw : float
    air_ext_lw : float
    air_ssa_lw : float
    veg_ssa_lw : float
    veg_fsd_const : float
    veg_contact_fraction_const : float
    ground_albedo_dir_mult_fact : float
    use_sw_direct_albedo : bool
    height : float array
    building_frac : float array
    veg_frac : float array
    building_scale : float array
    veg_scale : float array
    alb_roof : float array
    emis_roof : float array
    alb_wall : float array
    emis_wall : float array
    roof_albedo_dir_mult_fact : float array
    wall_specular_frac : float array
    radmeltfact : float
    raincover : float
    rainmaxres : float
    resp_a : float array
    resp_b : float array
    roughlenheatmethod : int
    roughlenmommethod : int
    runofftowater : float
    s1 : float
    s2 : float
    sathydraulicconduct : float array
    sddfull : float array
    sdd_id : float array
    smdmethod : int
    snowalb : float
    snowalbmax : float
    snowalbmin : float
    snowpacklimit : float array
    snowdens : float array
    snowdensmax : float
    snowdensmin : float
    snowfallcum : float
    snowfrac : float array
    snowlimbldg : float
    snowlimpaved : float
    snowpack : float array
    snowprof_24hr : float array
    snowuse : int
    soildepth : float array
    stabilitymethod : int
    startdls : int
    soilstore_surf : float array
    soilstorecap_surf : float array
    state_surf : float array
    statelimit_surf : float array
    wetthresh_surf : float array
    soilstore_roof : float array
    soilstorecap_roof : float array
    state_roof : float array
    statelimit_roof : float array
    wetthresh_roof : float array
    soilstore_wall : float array
    soilstorecap_wall : float array
    state_wall : float array
    statelimit_wall : float array
    wetthresh_wall : float array
    storageheatmethod : int
    storedrainprm : float array
    surfacearea : float
    tair_av : float
    tau_a : float
    tau_f : float
    tau_r : float
    baset_cooling : float array
    baset_heating : float array
    tempmeltfact : float
    th : float
    theta_bioco2 : float array
    timezone : float
    tl : float
    trafficrate : float array
    trafficunits : float
    sfr_surf : float array
    tsfc_roof : float array
    tsfc_wall : float array
    tsfc_surf : float array
    temp_roof : float array
    temp_wall : float array
    temp_surf : float array
    tin_roof : float array
    tin_wall : float array
    tin_surf : float array
    k_wall : float array
    k_roof : float array
    k_surf : float array
    cp_wall : float array
    cp_roof : float array
    cp_surf : float array
    dz_wall : float array
    dz_roof : float array
    dz_surf : float array
    tmin_id : float
    tmax_id : float
    lenday_id : float
    traffprof_24hr : float array
    ts5mindata_ir : float array
    tstep : int
    tstep_prev : int
    veg_type : int
    waterdist : float array
    waterusemethod : int
    wuday_id : float array
    decidcap_id : float
    albdectr_id : float
    albevetr_id : float
    albgrass_id : float
    porosity_id : float
    wuprofa_24hr : float array
    wuprofm_24hr : float array
    z : float
    z0m_in : float
    zdm_in : float
    dataoutblocksuews : float array
    dataoutblocksnow : float array
    dataoutblockestm : float array
    dataoutblockrsl : float array
    dataoutblockbeers : float array
    dataoutblockdebug : float array
    dataoutblockspartacus : float array
    dataoutblockestmext : float array
    dailystateblock : float array
    
    ================================================
     below is for debugging
     WRITE(year_txt, '(I4)') INT(iy)
     WRITE(id_text, '(I3)') INT(id)
     WRITE(it_text, '(I4)') INT(it)
     WRITE(imin_text, '(I4)') INT(imin)
     FileStateInit = './'//TRIM(ADJUSTL(year_txt))//'_'&
     //TRIM(ADJUSTL(id_text))//'_'&
     //TRIM(ADJUSTL(it_text))//'_'&
     //TRIM(ADJUSTL(imin_text))//'_'&
     //'state_init.nml'
     OPEN(12, file=FileStateInit, position='rewind')
     write(12, *) '&state_init'
     write(12, *) 'aerodynamicresistancemethod=', aerodynamicresistancemethod
     write(12, *) 'ah_min=', ah_min
     write(12, *) 'ahprof_24hr=', ahprof_24hr
     write(12, *) 'ah_slope_cooling=', ah_slope_cooling
     write(12, *) 'ah_slope_heating=', ah_slope_heating
     write(12, *) 'alb=', alb
     write(12, *) 'albmax_dectr=', albmax_dectr
     write(12, *) 'albmax_evetr=', albmax_evetr
     write(12, *) 'albmax_grass=', albmax_grass
     write(12, *) 'albmin_dectr=', albmin_dectr
     write(12, *) 'albmin_evetr=', albmin_evetr
     write(12, *) 'albmin_grass=', albmin_grass
     write(12, *) 'alpha_bioco2=', alpha_bioco2
     write(12, *) 'alpha_enh_bioco2=', alpha_enh_bioco2
     write(12, *) 'alt=', alt
     write(12, *) 'avkdn=', avkdn
     write(12, *) 'avrh=', avrh
     write(12, *) 'avu1=', avu1
     write(12, *) 'baset=', baset
     write(12, *) 'basete=', basete
     write(12, *) 'BaseT_HC=', BaseT_HC
     write(12, *) 'beta_bioco2=', beta_bioco2
     write(12, *) 'beta_enh_bioco2=', beta_enh_bioco2
     write(12, *) 'bldgh=', bldgh
     write(12, *) 'capmax_dec=', capmax_dec
     write(12, *) 'capmin_dec=', capmin_dec
     write(12, *) 'chanohm=', chanohm
     write(12, *) 'co2pointsource=', co2pointsource
     write(12, *) 'cpanohm=', cpanohm
     write(12, *) 'crwmax=', crwmax
     write(12, *) 'crwmin=', crwmin
     write(12, *) 'daywat=', daywat
     write(12, *) 'daywatper=', daywatper
     write(12, *) 'dectreeh=', dectreeh
     write(12, *) 'diagnose=', diagnose
     write(12, *) 'diagqn=', diagqn
     write(12, *) 'diagqs=', diagqs
     write(12, *) 'drainrt=', drainrt
     write(12, *) 'dt_since_start=', dt_since_start
     write(12, *) 'dqndt=', dqndt
     write(12, *) 'qn_av=', qn_av
     write(12, *) 'dqnsdt=', dqnsdt
     write(12, *) 'qn1_s_av=', qn1_s_av
     write(12, *) 'ef_umolco2perj=', ef_umolco2perj
     write(12, *) 'emis=', emis
     write(12, *) 'emissionsmethod=', emissionsmethod
     write(12, *) 'enef_v_jkm=', enef_v_jkm
     write(12, *) 'enddls=', enddls
     write(12, *) 'evetreeh=', evetreeh
     write(12, *) 'faibldg=', faibldg
     write(12, *) 'faidectree=', faidectree
     write(12, *) 'faievetree=', faievetree
     write(12, *) 'faut=', faut
     write(12, *) 'fcef_v_kgkm=', fcef_v_kgkm
     write(12, *) 'fcld_obs=', fcld_obs
     write(12, *) 'flowchange=', flowchange
     write(12, *) 'frfossilfuel_heat=', frfossilfuel_heat
     write(12, *) 'frfossilfuel_nonheat=', frfossilfuel_nonheat
     write(12, *) 'g1=', g1
     write(12, *) 'g2=', g2
     write(12, *) 'g3=', g3
     write(12, *) 'g4=', g4
     write(12, *) 'g5=', g5
     write(12, *) 'g6=', g6
     write(12, *) 'gdd_id=', gdd_id
     write(12, *) 'gddfull=', gddfull
     write(12, *) 'gridiv=', gridiv
     write(12, *) 'gsmodel=', gsmodel
     write(12, *) 'hdd_id=', hdd_id
     write(12, *) 'humactivity_24hr=', humactivity_24hr
     write(12, *) 'icefrac=', icefrac
     write(12, *) 'id=', id
     write(12, *) 'ie_a=', ie_a
     write(12, *) 'ie_end=', ie_end
     write(12, *) 'ie_m=', ie_m
     write(12, *) 'ie_start=', ie_start
     write(12, *) 'imin=', imin
     write(12, *) 'internalwateruse_h=', internalwateruse_h
     write(12, *) 'IrrFracEveTr=', IrrFracEveTr
     write(12, *) 'IrrFracDecTr=', IrrFracDecTr
     write(12, *) 'irrfracgrass=', irrfracgrass
     write(12, *) 'isec=', isec
     write(12, *) 'it=', it
     write(12, *) 'evapmethod=', evapmethod
     write(12, *) 'iy=', iy
     write(12, *) 'kkanohm=', kkanohm
     write(12, *) 'kmax=', kmax
     write(12, *) 'lai_id=', lai_id
     write(12, *) 'laicalcyes=', laicalcyes
     write(12, *) 'laimax=', laimax
     write(12, *) 'laimin=', laimin
     write(12, *) 'lai_obs=', lai_obs
     write(12, *) 'laipower=', laipower
     write(12, *) 'laitype=', laitype
     write(12, *) 'lat=', lat
     write(12, *) 'lenday_id=', lenday_id
     write(12, *) 'ldown_obs=', ldown_obs
     write(12, *) 'lng=', lng
     write(12, *) 'maxconductance=', maxconductance
     write(12, *) 'maxfcmetab=', maxfcmetab
     write(12, *) 'maxqfmetab=', maxqfmetab
     write(12, *) 'snowwater=', snowwater
     write(12, *) 'metforcingdata_grid=', metforcingdata_grid
     write(12, *) 'minfcmetab=', minfcmetab
     write(12, *) 'minqfmetab=', minqfmetab
     write(12, *) 'min_res_bioco2=', min_res_bioco2
     write(12, *) 'narp_emis_snow=', narp_emis_snow
     write(12, *) 'narp_trans_site=', narp_trans_site
     write(12, *) 'netradiationmethod=', netradiationmethod
     write(12, *) 'ohm_coef=', ohm_coef
     write(12, *) 'ohmincqf=', ohmincqf
     write(12, *) 'ohm_threshsw=', ohm_threshsw
     write(12, *) 'ohm_threshwd=', ohm_threshwd
     write(12, *) 'pipecapacity=', pipecapacity
     write(12, *) 'popdensdaytime=', popdensdaytime
     write(12, *) 'popdensnighttime=', popdensnighttime
     write(12, *) 'popprof_24hr=', popprof_24hr
     write(12, *) 'pormax_dec=', pormax_dec
     write(12, *) 'pormin_dec=', pormin_dec
     write(12, *) 'precip=', precip
     write(12, *) 'preciplimit=', preciplimit
     write(12, *) 'preciplimitalb=', preciplimitalb
     write(12, *) 'press_hpa=', press_hpa
     write(12, *) 'qf0_beu=', qf0_beu
     write(12, *) 'qf_a=', qf_a
     write(12, *) 'qf_b=', qf_b
     write(12, *) 'qf_c=', qf_c
     write(12, *) 'qn1_obs=', qn1_obs
     write(12, *) 'qh_obs=', qh_obs
     write(12, *) 'qs_obs=', qs_obs
     write(12, *) 'qf_obs=', qf_obs
     write(12, *) 'radmeltfact=', radmeltfact
     write(12, *) 'raincover=', raincover
     write(12, *) 'rainmaxres=', rainmaxres
     write(12, *) 'resp_a=', resp_a
     write(12, *) 'resp_b=', resp_b
     write(12, *) 'roughlenheatmethod=', roughlenheatmethod
     write(12, *) 'roughlenmommethod=', roughlenmommethod
     write(12, *) 'runofftowater=', runofftowater
     write(12, *) 's1=', s1
     write(12, *) 's2=', s2
     write(12, *) 'sathydraulicconduct=', sathydraulicconduct
     write(12, *) 'sddfull=', sddfull
     write(12, *) 'sdd_id=', sdd_id
     write(12, *) 'sfr_surf=', sfr_surf
     write(12, *) 'smdmethod=', smdmethod
     write(12, *) 'snowalb=', snowalb
     write(12, *) 'snowalbmax=', snowalbmax
     write(12, *) 'snowalbmin=', snowalbmin
     write(12, *) 'snowpacklimit=', snowpacklimit
     write(12, *) 'snowdens=', snowdens
     write(12, *) 'snowdensmax=', snowdensmax
     write(12, *) 'snowdensmin=', snowdensmin
     write(12, *) 'snowfallcum=', snowfallcum
     write(12, *) 'snowfrac=', snowfrac
     write(12, *) 'snowlimbldg=', snowlimbldg
     write(12, *) 'snowlimpaved=', snowlimpaved
     write(12, *) 'snowfrac_obs=', snowfrac_obs
     write(12, *) 'snowpack=', snowpack
     write(12, *) 'snowprof_24hr=', snowprof_24hr
     write(12, *) 'SnowUse=', SnowUse
     write(12, *) 'soildepth=', soildepth
     write(12, *) 'soilstore_id=', soilstore_id
     write(12, *) 'soilstorecap=', soilstorecap
     write(12, *) 'stabilitymethod=', stabilitymethod
     write(12, *) 'startdls=', startdls
     write(12, *) 'state_id=', state_id
     write(12, *) 'statelimit=', statelimit
     write(12, *) 'storageheatmethod=', storageheatmethod
     write(12, *) 'storedrainprm=', storedrainprm
     write(12, *) 'surfacearea=', surfacearea
     write(12, *) 'tair_av=', tair_av
     write(12, *) 'tau_a=', tau_a
     write(12, *) 'tau_f=', tau_f
     write(12, *) 'tau_r=', tau_r
     write(12, *) 'tmax_id=', tmax_id
     write(12, *) 'tmin_id=', tmin_id
     write(12, *) 'BaseT_Cooling=', BaseT_Cooling
     write(12, *) 'BaseT_Heating=', BaseT_Heating
     write(12, *) 'temp_c=', temp_c
     write(12, *) 'tempmeltfact=', tempmeltfact
     write(12, *) 'th=', th
     write(12, *) 'theta_bioco2=', theta_bioco2
     write(12, *) 'timezone=', timezone
     write(12, *) 'tl=', tl
     write(12, *) 'trafficrate=', trafficrate
     write(12, *) 'trafficunits=', trafficunits
     write(12, *) 'traffprof_24hr=', traffprof_24hr
     write(12, *) 'ts5mindata_ir=', ts5mindata_ir
     write(12, *) 'tstep=', tstep
     write(12, *) 'tstep_prev=', tstep_prev
     write(12, *) 'veg_type=', veg_type
     write(12, *) 'waterdist=', waterdist
     write(12, *) 'waterusemethod=', waterusemethod
     write(12, *) 'wetthresh=', wetthresh
     write(12, *) 'wu_m3=', wu_m3
     write(12, *) 'wuday_id=', wuday_id
     write(12, *) 'decidcap_id=', decidcap_id
     write(12, *) 'albdectr_id=', albdectr_id
     write(12, *) 'albevetr_id=', albevetr_id
     write(12, *) 'albgrass_id=', albgrass_id
     write(12, *) 'porosity_id=', porosity_id
     write(12, *) 'wuprofa_24hr=', wuprofa_24hr
     write(12, *) 'wuprofm_24hr=', wuprofm_24hr
     write(12, *) 'xsmd=', xsmd
     write(12, *) 'z=', z
     write(12, *) 'z0m_in=', z0m_in
     write(12, *) 'zdm_in=', zdm_in
     write(12, *) '/'
     WRITE(12, *) ''
     CLOSE(12)
    ================================================
    """
    _suews_driver.f90wrap_suews_cal_multitsteps(metforcingblock=metforcingblock, \
        len_sim=len_sim, ah_min=ah_min, ahprof_24hr=ahprof_24hr, \
        ah_slope_cooling=ah_slope_cooling, ah_slope_heating=ah_slope_heating, \
        alb=alb, albmax_dectr=albmax_dectr, albmax_evetr=albmax_evetr, \
        albmax_grass=albmax_grass, albmin_dectr=albmin_dectr, \
        albmin_evetr=albmin_evetr, albmin_grass=albmin_grass, \
        alpha_bioco2=alpha_bioco2, alpha_enh_bioco2=alpha_enh_bioco2, alt=alt, \
        baset=baset, basete=basete, beta_bioco2=beta_bioco2, \
        beta_enh_bioco2=beta_enh_bioco2, bldgh=bldgh, capmax_dec=capmax_dec, \
        capmin_dec=capmin_dec, chanohm=chanohm, co2pointsource=co2pointsource, \
        cpanohm=cpanohm, crwmax=crwmax, crwmin=crwmin, daywat=daywat, \
        daywatper=daywatper, dectreeh=dectreeh, diagmethod=diagmethod, \
        diagnose=diagnose, drainrt=drainrt, dt_since_start=dt_since_start, \
        dqndt=dqndt, qn_av=qn_av, dqnsdt=dqnsdt, qn_s_av=qn_s_av, \
        ef_umolco2perj=ef_umolco2perj, emis=emis, emissionsmethod=emissionsmethod, \
        enef_v_jkm=enef_v_jkm, enddls=enddls, evetreeh=evetreeh, faibldg=faibldg, \
        faidectree=faidectree, faievetree=faievetree, faut=faut, \
        fcef_v_kgkm=fcef_v_kgkm, flowchange=flowchange, \
        frfossilfuel_heat=frfossilfuel_heat, \
        frfossilfuel_nonheat=frfossilfuel_nonheat, g_max=g_max, g_k=g_k, \
        g_q_base=g_q_base, g_q_shape=g_q_shape, g_t=g_t, g_sm=g_sm, gdd_id=gdd_id, \
        gddfull=gddfull, gridiv=gridiv, gsmodel=gsmodel, h_maintain=h_maintain, \
        hdd_id=hdd_id, humactivity_24hr=humactivity_24hr, icefrac=icefrac, \
        ie_a=ie_a, ie_end=ie_end, ie_m=ie_m, ie_start=ie_start, \
        internalwateruse_h=internalwateruse_h, irrfracpaved=irrfracpaved, \
        irrfracbldgs=irrfracbldgs, irrfracevetr=irrfracevetr, \
        irrfracdectr=irrfracdectr, irrfracgrass=irrfracgrass, \
        irrfracbsoil=irrfracbsoil, irrfracwater=irrfracwater, kkanohm=kkanohm, \
        kmax=kmax, lai_id=lai_id, laimax=laimax, laimin=laimin, laipower=laipower, \
        laitype=laitype, lat=lat, lng=lng, maxconductance=maxconductance, \
        maxfcmetab=maxfcmetab, maxqfmetab=maxqfmetab, snowwater=snowwater, \
        minfcmetab=minfcmetab, minqfmetab=minqfmetab, min_res_bioco2=min_res_bioco2, \
        narp_emis_snow=narp_emis_snow, narp_trans_site=narp_trans_site, \
        netradiationmethod=netradiationmethod, ohm_coef=ohm_coef, ohmincqf=ohmincqf, \
        ohm_threshsw=ohm_threshsw, ohm_threshwd=ohm_threshwd, \
        pipecapacity=pipecapacity, popdensdaytime=popdensdaytime, \
        popdensnighttime=popdensnighttime, popprof_24hr=popprof_24hr, \
        pormax_dec=pormax_dec, pormin_dec=pormin_dec, preciplimit=preciplimit, \
        preciplimitalb=preciplimitalb, qf0_beu=qf0_beu, qf_a=qf_a, qf_b=qf_b, \
        qf_c=qf_c, nlayer=nlayer, \
        n_vegetation_region_urban=n_vegetation_region_urban, \
        n_stream_sw_urban=n_stream_sw_urban, n_stream_lw_urban=n_stream_lw_urban, \
        sw_dn_direct_frac=sw_dn_direct_frac, air_ext_sw=air_ext_sw, \
        air_ssa_sw=air_ssa_sw, veg_ssa_sw=veg_ssa_sw, air_ext_lw=air_ext_lw, \
        air_ssa_lw=air_ssa_lw, veg_ssa_lw=veg_ssa_lw, veg_fsd_const=veg_fsd_const, \
        veg_contact_fraction_const=veg_contact_fraction_const, \
        ground_albedo_dir_mult_fact=ground_albedo_dir_mult_fact, \
        use_sw_direct_albedo=use_sw_direct_albedo, height=height, \
        building_frac=building_frac, veg_frac=veg_frac, \
        building_scale=building_scale, veg_scale=veg_scale, alb_roof=alb_roof, \
        emis_roof=emis_roof, alb_wall=alb_wall, emis_wall=emis_wall, \
        roof_albedo_dir_mult_fact=roof_albedo_dir_mult_fact, \
        wall_specular_frac=wall_specular_frac, radmeltfact=radmeltfact, \
        raincover=raincover, rainmaxres=rainmaxres, resp_a=resp_a, resp_b=resp_b, \
        roughlenheatmethod=roughlenheatmethod, roughlenmommethod=roughlenmommethod, \
        runofftowater=runofftowater, s1=s1, s2=s2, \
        sathydraulicconduct=sathydraulicconduct, sddfull=sddfull, sdd_id=sdd_id, \
        smdmethod=smdmethod, snowalb=snowalb, snowalbmax=snowalbmax, \
        snowalbmin=snowalbmin, snowpacklimit=snowpacklimit, snowdens=snowdens, \
        snowdensmax=snowdensmax, snowdensmin=snowdensmin, snowfallcum=snowfallcum, \
        snowfrac=snowfrac, snowlimbldg=snowlimbldg, snowlimpaved=snowlimpaved, \
        snowpack=snowpack, snowprof_24hr=snowprof_24hr, snowuse=snowuse, \
        soildepth=soildepth, stabilitymethod=stabilitymethod, startdls=startdls, \
        soilstore_surf=soilstore_surf, soilstorecap_surf=soilstorecap_surf, \
        state_surf=state_surf, statelimit_surf=statelimit_surf, \
        wetthresh_surf=wetthresh_surf, soilstore_roof=soilstore_roof, \
        soilstorecap_roof=soilstorecap_roof, state_roof=state_roof, \
        statelimit_roof=statelimit_roof, wetthresh_roof=wetthresh_roof, \
        soilstore_wall=soilstore_wall, soilstorecap_wall=soilstorecap_wall, \
        state_wall=state_wall, statelimit_wall=statelimit_wall, \
        wetthresh_wall=wetthresh_wall, storageheatmethod=storageheatmethod, \
        storedrainprm=storedrainprm, surfacearea=surfacearea, tair_av=tair_av, \
        tau_a=tau_a, tau_f=tau_f, tau_r=tau_r, baset_cooling=baset_cooling, \
        baset_heating=baset_heating, tempmeltfact=tempmeltfact, th=th, \
        theta_bioco2=theta_bioco2, timezone=timezone, tl=tl, \
        trafficrate=trafficrate, trafficunits=trafficunits, sfr_surf=sfr_surf, \
        tsfc_roof=tsfc_roof, tsfc_wall=tsfc_wall, tsfc_surf=tsfc_surf, \
        temp_roof=temp_roof, temp_wall=temp_wall, temp_surf=temp_surf, \
        tin_roof=tin_roof, tin_wall=tin_wall, tin_surf=tin_surf, k_wall=k_wall, \
        k_roof=k_roof, k_surf=k_surf, cp_wall=cp_wall, cp_roof=cp_roof, \
        cp_surf=cp_surf, dz_wall=dz_wall, dz_roof=dz_roof, dz_surf=dz_surf, \
        tmin_id=tmin_id, tmax_id=tmax_id, lenday_id=lenday_id, \
        traffprof_24hr=traffprof_24hr, ts5mindata_ir=ts5mindata_ir, tstep=tstep, \
        tstep_prev=tstep_prev, veg_type=veg_type, waterdist=waterdist, \
        waterusemethod=waterusemethod, wuday_id=wuday_id, decidcap_id=decidcap_id, \
        albdectr_id=albdectr_id, albevetr_id=albevetr_id, albgrass_id=albgrass_id, \
        porosity_id=porosity_id, wuprofa_24hr=wuprofa_24hr, \
        wuprofm_24hr=wuprofm_24hr, z=z, z0m_in=z0m_in, zdm_in=zdm_in, \
        dataoutblocksuews=dataoutblocksuews, dataoutblocksnow=dataoutblocksnow, \
        dataoutblockestm=dataoutblockestm, dataoutblockrsl=dataoutblockrsl, \
        dataoutblockbeers=dataoutblockbeers, dataoutblockdebug=dataoutblockdebug, \
        dataoutblockspartacus=dataoutblockspartacus, \
        dataoutblockestmext=dataoutblockestmext, dailystateblock=dailystateblock)

def suews_cal_sunposition(year, idectime, utc, locationlatitude, \
    locationlongitude, locationaltitude):
    """
    sunazimuth, sunzenith = suews_cal_sunposition(year, idectime, utc, \
        locationlatitude, locationlongitude, locationaltitude)
    
    
    Defined at suews_ctrl_driver.fpp lines 4413-4420
    
    Parameters
    ----------
    year : float
    idectime : float
    utc : float
    locationlatitude : float
    locationlongitude : float
    locationaltitude : float
    
    Returns
    -------
    sunazimuth : float
    sunzenith : float
    
    """
    sunazimuth, sunzenith = _suews_driver.f90wrap_suews_cal_sunposition(year=year, \
        idectime=idectime, utc=utc, locationlatitude=locationlatitude, \
        locationlongitude=locationlongitude, locationaltitude=locationaltitude)
    return sunazimuth, sunzenith

def cal_tair_av(tair_av_prev, dt_since_start, tstep, temp_c):
    """
    tair_av_next = cal_tair_av(tair_av_prev, dt_since_start, tstep, temp_c)
    
    
    Defined at suews_ctrl_driver.fpp lines 4427-4448
    
    Parameters
    ----------
    tair_av_prev : float
    dt_since_start : int
    tstep : int
    temp_c : float
    
    Returns
    -------
    tair_av_next : float
    
    """
    tair_av_next = _suews_driver.f90wrap_cal_tair_av(tair_av_prev=tair_av_prev, \
        dt_since_start=dt_since_start, tstep=tstep, temp_c=temp_c)
    return tair_av_next

def cal_tsfc(qh, avdens, avcp, ra, temp_c):
    """
    tsfc_c = cal_tsfc(qh, avdens, avcp, ra, temp_c)
    
    
    Defined at suews_ctrl_driver.fpp lines 4450-4460
    
    Parameters
    ----------
    qh : float
    avdens : float
    avcp : float
    ra : float
    temp_c : float
    
    Returns
    -------
    tsfc_c : float
    
    """
    tsfc_c = _suews_driver.f90wrap_cal_tsfc(qh=qh, avdens=avdens, avcp=avcp, ra=ra, \
        temp_c=temp_c)
    return tsfc_c


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "suews_driver".')

for func in _dt_array_initialisers:
    func()
