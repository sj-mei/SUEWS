!
!
!
!
  module modulestebbsprecision
!
  use iso_fortran_env, only : real64   
!
  implicit none
!
  integer,parameter :: rprc = real64
!
  end module
!
!
!
!
  module modulestebbs
!
  use modulestebbsprecision
!
  real(rprc),parameter :: sigma = 5.670E-8
!
  integer,save :: flgtimecheck = 1
  integer :: resolution
  integer :: time_st, time_ed, count_p_sec, count_max ! Time check
!
  integer :: nbtype
!
  character(len=256),allocatable,dimension(:) :: fnmls, cases
!
  type :: LBM
!
    character(len=256) ::                          &
    BuildingType,                                  &
    BuildingName,                                  &
    fnmlLBM,                                       &
    case
!
    integer         :: idLBM
    integer         :: flginit = 0
    integer         :: appliance_totalnumber
!
    real(rprc) ::                                  &
    Qtotal_heating,                                &
    Qtotal_cooling,                                &
    Qmetabolic_sensible,                           &
    Qmetabolic_latent,                             &
    Qtotal_water_tank,                             &
    qhwtDrain,                                     &
    ratio_window_wall,                             &
    Afootprint,                                    &
    height_building,                               &
    wallExternalArea,                              &
    ratioInternalVolume,                           &
    thickness_wallroof,                            &
    thickness_groundfloor,                         &
    depth_ground,                                  &
    thickness_window,                              &
    conv_coeff_intwallroof,                        &
    conv_coeff_indoormass,                         &
    conv_coeff_intgroundfloor,                     &
    conv_coeff_intwindow,                          &
    conv_coeff_extwallroof,                        &
    conv_coeff_extwindow,                          &
    conductivity_wallroof,                         &
    conductivity_groundfloor,                      &
    conductivity_window,                           &
    conductivity_ground,                           &
    density_wallroof,                              &
    weighting_factor_heatcapacity_wallroof,        &
    density_groundfloor,                           &
    density_window,                                &
    density_indoormass,                            &
    density_air_ind,                               &
    cp_wallroof,                                   &
    cp_groundfloor,                                &
    cp_window,                                     &
    cp_indoormass,                                 &
    cp_air_ind,                                    &
    emissivity_extwallroof,                        &
    emissivity_intwallroof,                        &
    emissivity_indoormass,                         &
    emissivity_extwindow,                          &
    emissivity_intwindow,                          &
    windowTransmissivity,                          &
    windowAbsorbtivity,                            &
    windowReflectivity,                            &
    wallTransmisivity,                             &
    wallAbsorbtivity,                              &
    wallReflectivity,                              &
    BVF_extwall,                                   &
    GVF_extwall,                                   &
    SVF_extwall,                                   &
    occupants,                                     &
    metabolic_rate,                                &
    ratio_metabolic_latent_sensible,               &
    appliance_power_rating,                        &
    appliance_usage_factor,                        &
    maxheatingpower_air,                           &
    heating_efficiency_air,                        &
    maxcoolingpower_air,                           &
    coeff_performance_cooling,                     &
    Vair_ind,                                      &
    ventilation_rate,                              &
    Awallroof,                                     &
    Vwallroof,                                     &
    Vgroundfloor,                                  &
    Awindow,                                       &
    Vwindow,                                       &
    Vindoormass,                                   &
    Aindoormass,                                   &
    Tair_ind,                                      &
    Tindoormass,                                   &
    Tintwallroof,                                  &
    Textwallroof,                                  &
    Tintwindow,                                    &
    Textwindow,                                    &
    Tintgroundfloor,                               &
    Textgroundfloor,                               &
    Twater_tank,                                   &
    Tintwall_tank,                                 &
    Textwall_tank,                                 &
    thickness_tankwall,                            &
    Tincomingwater_tank,                           &
    Vwater_tank,                                   &
    Asurf_tank,                                    &
    Vwall_tank,                                    &
    setTwater_tank,                                &
    init_wtTs,                                     &
    Twater_vessel,                                 &
    Tintwall_vessel,                               &
    Textwall_vessel,                               &
    thickness_wall_vessel,                         &
    Vwater_vessel,                                 &
    Awater_vessel,                                 &
    Vwall_vessel,                                  &
    flowrate_water_supply,                         &
    flowrate_water_drain,                          &
    single_flowrate_water_supply,                  &
    single_flowrate_water_drain,                   &
    cp_water,                                      &
    cp_water_tank,                                 &
    cp_wall_tank,                                  &
    cp_wall_vessel,                                &
    density_water,                                 &
    density_wall_tank,                             &
    density_wall_vessel,                           &
    BVF_tank,                                      &
    MVF_tank,                                      &
    conductivity_wall_tank,                        &
    conv_coeff_intwall_tank,                       &
    conv_coeff_extwall_tank,                       &
    emissivity_extwall_tank,                       &
    conductivity_wall_vessel,                      &
    conv_coeff_intwall_vessel,                     &
    conv_coeff_extwall_vessel,                     &
    emissivity_extwall_vessel,                     &
    maxheatingpower_water,                         &
    heating_efficiency_water,                      &
    minVwater_vessel,                              &
    minHeatingPower_DHW,                           &
    HeatingPower_DHW
!
!
!
!   STEBBS output
!
    real(rprc) ::       &
    qfm_dom,            & ! Metabolic sensible and latent heat
    qheat_dom,          & ! Hourly heating load  [W]
    qcool_dom,          & ! Hourly cooling load  [W]
    qfb_hw_dom,         & ! Hot water
    qfb_dom_air,        & ! Sensible heat to air [W]
    dom_temp,           & ! Domain temperature   [W]
    QStar,              & ! Net radiation        [W m-2]
    QEC,                & ! Energy use           [W m-2]
    QH,                 & ! Sensible heat flux   [W m-2]
    QS,                 & ! Storage heat flux    [W m-2]
    QBAE,               & ! Building exchange    [W m-2]
    QWaste                ! Waste heating        [W m-2]
!
!
!
    real(rprc),dimension(2) :: Ts, initTs
    real(rprc),dimension(4) :: h_i, k_eff
    real(rprc),dimension(2) :: h_o
    real(rprc),dimension(5) :: rho
    real(rprc),dimension(5) :: Cp
    real(rprc),dimension(5) :: emis
    real(rprc),dimension(3) :: wiTAR, waTAR
    real(rprc),dimension(3) :: viewFactors
    real(rprc),dimension(3) :: occupantData 
    real(rprc),dimension(3) :: HTsAverage, HWTsAverage
    real(rprc),dimension(3) :: HWPowerAverage
!
    real(rprc),dimension(25) :: EnergyExchanges = 0.0
!
!
!
  end type
!
  type(LBM),allocatable,dimension(:) :: blds
!
  end module modulestebbs
!
!
!
!
  module modulestebbsfunc 
!
  use modulestebbsprecision
!
  implicit none
!
  contains
!
!
!
!/**************** Water lost to drains *****************/
!// rho - density of water
!// Cp - specific heat capacity
!// vFR - volume Flow Rate
!// Tout - temperature (K) of water in vessel lost to drains
  function waterUseEnergyLossToDrains(rho,Cp,vFRo,Tout,timeResolution) result (q_wt)
!
  use modulestebbsprecision
!
  implicit none
!
  integer,intent(in)          :: timeResolution
  real(rprc),intent(in) :: rho, Cp, vFRo, Tout
  real(rprc)            :: q_wt
!
  q_wt = rho*Cp*Tout*(vFRo*timeResolution) ! // kg/m3 . J/K.kg . K . (m3/s . s) (units cancel to J)
!
  end function
!
!
!
!
  function indoorConvectionHeatTransfer(h,A,Twi,Ti) result (ind_cht)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: h, A, Twi, Ti
  real(rprc)             :: ind_cht
!
  ind_cht = h * A * (Ti - Twi)
!
  end function indoorConvectionHeatTransfer
!
!
!
!
  function internalConvectionHeatTransfer(h,A,Tio,Ti) result (int_cht)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: h, A, Tio, Ti
  real(rprc)             :: int_cht
!
  int_cht = h * A * (Ti - Tio)
!
  end function internalConvectionHeatTransfer
!
!
!
!
  function indoorRadiativeHeatTransfer() result (q)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc) :: q
!
  q = 0.0
!
  return
!
  end function indoorRadiativeHeatTransfer
!
!
!
!
  function outdoorConvectionHeatTransfer(h,A,Two,Ta) result (out_cht)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: h, A, Two, Ta
  real(rprc)             :: out_cht
!
  out_cht = h * A * (Two - Ta)
!
  return
!
  end function outdoorConvectionHeatTransfer
!
!
!
!
  function outdoorRadiativeHeatTransfer(f,A,emis,Two,Ts) result (q)
!
  use modulestebbsprecision
  use modulestebbs, only : sigma
!
  implicit none
!
  real(rprc),intent(in)  :: f,A,emis,Two,Ts
  real(rprc)             :: q
!
  q = A * f * sigma * emis * (Two ** 4.0 - Ts ** 4.0 )
!
  return
!
  end function outdoorRadiativeHeatTransfer
!
!
!
!
  function lwoutdoorRadiativeHeatTransfer(A,emis,Two,lw) result (q)
!
  use modulestebbsprecision
  use modulestebbs, only : sigma
!
  implicit none
!
  real(rprc),intent(in)  :: A,emis,Two,lw
  real(rprc)             :: q
!
  q = A * sigma * emis * (Two ** 4.0) - emis*lw*a ! Revised based on Yiqing's discovery
!
  return
!
  end function lwoutdoorRadiativeHeatTransfer
!
!
!
!/************************************************************/
!// Window Solar Insolation
!// Kwall - Irradiance incident on vertical window/wall
!// Tr - Effective Transmissivity
!// NOTE: This should be added as heat gain to internal single mass object
!
  function windowInsolation(Irr,Tr,A) result(wi_in)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: Irr,Tr,A
  real(rprc)             :: wi_in
!
  wi_in = Irr * Tr * A
!
  return
!
  end function windowInsolation
!
!
!
!/************************************************************/
!//  Solar Insolation on surface
!// Irr - Irradiance incident horiozntal roof or vertical wall (irr2)
!// Ab - Effective Absorptance of surface
!// NOTE: This should be added as heat gain to external surface of wall
!
  function wallInsolation(Irr,Ab,A) result (wa_in)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: Irr,Ab,A
  real(rprc)             :: wa_in
!
  wa_in = Irr * Ab * A
!
  return
!
  end function wallInsolation
!
!
!
!
  function wallConduction(k_eff,A,Twi,Two,L) result (wa_co)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: k_eff,Twi,Two,A,L
  real(rprc)             :: wa_co 
!
  wa_co = k_eff * A *((Twi - Two)/L)
!
  end function wallConduction
!
!
!
!
  function windowConduction(k_eff,A,Twi,Two,L) result(wi_co)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: k_eff,Twi,Two,A,L
  real(rprc)             :: wi_co
!
  wi_co = k_eff * A *((Twi - Two)/L)
!
  end function windowConduction
!
!
!
!
  function heating(Ts,Ti,epsilon,P) result (q_heating)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: Ts,Ti,epsilon,P
  real(rprc)             :: q_heating
!
  q_heating = 0.0
!
  if (Ti < Ts ) then
      q_heating = (P - (P/exp(Ts-Ti)))*epsilon
  endif
!
  end function heating
!
!
!
!
  function ventilationHeatTransfer(rho,Cp,V,To,Ti) result(q_in)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: rho,Cp,V,To,Ti
  real(rprc)             :: q_in
!
  q_in = rho*Cp*V*(To-Ti)
! 
  end function ventilationHeatTransfer
!
!
!
!
  function additionalSystemHeatingEnergy(q_heating,epsilon) result(qH_additional)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: q_heating,epsilon
  real(rprc)             :: qH_additional
!
  qH_additional = 0.0
  qH_additional = (q_heating/epsilon) - q_heating
! 
  end function additionalSystemHeatingEnergy
!
!
!
!
  function cooling(Ts,Ti,COP,P) result(q_cooling)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: Ts,Ti,COP,P
  real(rprc)             :: q_cooling
!
  q_cooling = 0.0
!
  if (Ti > Ts) then
      q_cooling = P - (P/exp(Ti-Ts))
  endif
!
  end function cooling
!
!
!
!
  function additionalSystemCoolingEnergy(q_cooling,COP) result(qC_additional)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)  :: q_cooling,COP
  real(rprc)             :: qC_additional
!
  qC_additional = q_cooling/COP
!
  end function additionalSystemCoolingEnergy
!
!
!
!
  function internalOccupancyGains(Occupants,metRate,LSR) result(qSL)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in)   :: Occupants,metRate,LSR
  real(rprc)              :: qSen,qLat
  real(rprc),dimension(2) :: qSL
!
!  qSen = (metRate*Occupants)/(1.0+LSR)
!  qLat = (metRate*Occupants)*LSR/(1.0+LSR)
!
  qSL(1) = (metRate*Occupants)/(1.0+LSR) 
  qSL(2) = (metRate*Occupants)*LSR/(1.0+LSR)
!
  end function internalOccupancyGains
!
!
!
!
  function internalApplianceGains(P,f,n) result (qapp)
!
  use modulestebbsprecision
!
  implicit none
!
  integer,intent(in)           :: n
  real(rprc),intent(in)  :: P,f
  real(rprc)             :: qapp
!
  qapp = P*f*n
!
  end function internalApplianceGains
!
!
!
!
  function ext_conv_coeff(wind_speed, dT) result (hc)
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc),intent(in) :: wind_speed, dT
  real(rprc)            :: hn, a, b, Rf, hcglass, hc
!
    hn = 1.31 * (abs(dT) ** (1.0/3.0))
    a = 3.26
    b = 0.89
    Rf = 1.67 ! # Rough brick used from E+ Engineering Reference guide, Walton 1981
    hcglass = ((hn ** 2)+((a*(wind_speed ** b))** 2))** (0.5)
    hc = hn + Rf*(hcglass - hn)
!
  end function ext_conv_coeff
!
!
!
  end module modulestebbsfunc
!
!
!
!
  module modulesuewsstebbscouple
!
  use modulestebbsprecision
!
  implicit none
!
  real(rprc) ::       Tair_out, Tsurf, Tground_deep,  &
                      density_air_out, cp_air_out,    &
                      Qsw_dn_extroof, Qsw_dn_extwall, &
                      Qlw_dn_extwall, Qlw_dn_extroof
!
!
!
  type :: suewsprop
!
    integer                                     :: ntstep, timestep
    character(len=256),allocatable,dimension(:) :: datetime, hourmin
    real(rprc),allocatable,dimension(:)   :: Tair,Tsurf,Kwall,Kroof,ws,Lroof,Lwall
!
    integer                                     :: ntskip
    character(len=256),allocatable,dimension(:) :: datetime_exch, hourmin_exch
    real(rprc),allocatable,dimension(:)   :: Tair_exch,Tsurf_exch,Kwall_exch,Kroof_exch,ws_exch,Lroof_exch,Lwall_exch
!
  end type
!
!
!
  type(suewsprop) :: sout
!
  end module modulesuewsstebbscouple
!
!
!
!
  subroutine setdatetime(datetimeLine)
!
  use modulestebbsprecision
  use modulesuewsstebbscouple,only : sout
!
  implicit none
!
  real(rprc),dimension(5),intent(in) :: datetimeLine
!
  integer          :: i
  character(len=4) :: cyear
  character(len=2) :: cmonth, cday, chour, cmin, csec
  integer,dimension(12) :: stmonth
  integer,dimension(12) :: stmonth_nonleap = (/0,31,59,90,120,151,181,212,243,273,304,334/)
  integer,dimension(12) :: stmonth_leap = (/0,31,60,91,121,152,182,213,244,274,305,335/)
!
!
! Year
!
  write(cyear,'(i4)')int(datetimeLine(1))
!
  if ( mod(int(datetimeLine(1)),4) == 0 ) then
      stmonth = stmonth_leap
  else
      stmonth = stmonth_nonleap
  endif
!
! DOY
!
  do i = 1,11,1 
      if (stmonth(i) < datetimeLine(2) .and. datetimeLine(2) <= stmonth(i+1) ) then
          write(cmonth,'(i2.2)')i
          write(cday,'(i2.2)')int(datetimeLine(2)) - stmonth(i)
      endif
      if (stmonth(12) < datetimeLine(2) ) then
          write(cmonth,'(i2.2)')12
          write(cday,'(i2.2)')int(datetimeLine(2)) - stmonth(12)
      endif
  enddo
!
!
!
  write(chour,'(i2.2)') int(datetimeLine(3))
  write(cmin,'(i2.2)') int(datetimeLine(4))
  write(csec,'(i2.2)') 0
!
  sout%datetime(1) = trim(cyear//'-'//cmonth//'-'//cday)
  sout%hourmin(1) = trim(chour//':'//cmin//':'//csec)
!
!
!
  return
!
  end subroutine setdatetime
!
!
!
!
  subroutine stebbsonlinecouple(timestep, datetimeLine, Tair_sout,Tsurf_sout, &
                          Kroof_sout,Kwall_sout,Lwall_sout,Lroof_sout,ws)
!
  use modulestebbs
  use modulesuewsstebbscouple
!
  implicit none
!
  integer      :: i
  integer,intent(in) :: timestep
  integer,save :: flginit = 0
!
  real(rprc),intent(in) :: Tair_sout,Tsurf_sout,Kroof_sout, &
                           Kwall_sout,Lwall_sout,Lroof_sout,ws
  real(rprc),dimension(5),intent(in) :: datetimeLine
!
  namelist/settings/nbtype,resolution
  namelist/io/cases
!
!
!
  if ( flginit == 0 ) then
!
      open(8,file='./RunControl_STEBBS.nml',status='old',form='formatted')
      read(8,nml=settings)
      allocate(cases(nbtype))
      read(8,nml=io)
      close(8)
!
      allocate(fnmls(nbtype))
      allocate(blds(nbtype))
!
      write(*,*) '++++ SUEWS-STEBBS coupling'
      write(*,*) '    + Total building type : ' , nbtype
      do i = 1, nbtype, 1
          write(*,*) '    + Cases title         : ' , i , trim(cases(i))
      enddo
!
      if ( resolution <= 0 ) resolution = timestep
!
!
!
      do i = 1, nbtype, 1
          fnmls(i) = './BuildClasses/' // trim(cases(i)) // '.nml'
          call create_building(cases(i),blds(i),i)
      enddo
!
!
!
      sout%ntstep = 1
      allocate(sout%datetime(sout%ntstep))
      allocate(sout%hourmin(sout%ntstep))
      allocate(sout%Tair(   sout%ntstep))
      allocate(sout%Tsurf(  sout%ntstep))
      allocate(sout%Kwall(  sout%ntstep))
      allocate(sout%Kroof(  sout%ntstep))
      allocate(sout%ws(     sout%ntstep))
      allocate(sout%Lroof(  sout%ntstep))
      allocate(sout%Lwall(  sout%ntstep))
      allocate(sout%datetime_exch(sout%ntstep))
      allocate(sout%hourmin_exch(sout%ntstep))
      allocate(sout%Tair_exch(   sout%ntstep))
      allocate(sout%Tsurf_exch(  sout%ntstep))
      allocate(sout%Kwall_exch(  sout%ntstep))
      allocate(sout%Kroof_exch(  sout%ntstep))
      allocate(sout%ws_exch(     sout%ntstep))
      allocate(sout%Lroof_exch(  sout%ntstep))
      allocate(sout%Lwall_exch(  sout%ntstep))
!
  endif
!
!
!
! Hand over SUEWS output to STEBBS input
!
  sout%Tair(1)  = Tair_sout 
  sout%Tsurf(1) = Tsurf_sout
  sout%Kroof(1) = Kroof_sout
  sout%Kwall(1) = Kwall_sout
  sout%Lwall(1) = Lwall_sout
  sout%Lroof(1) = Lroof_sout
  sout%timestep = timestep
  sout%Tair_exch(1)  = Tair_sout
  sout%Tsurf_exch(1) = Tsurf_sout
  sout%ws_exch(1)    = ws
!
!
!
  call setdatetime(datetimeLine)
!
!
! Time integration for each building type
!
  do i = 1, nbtype, 1
      call suewsstebbscouple(blds(i))
  enddo
!
!
!
! Mush-up building-wise output to cast back to SUEWS
!
! SHOULD DO THIS HERE
!
!
!
  flginit = 1
!
  return
!
  end subroutine stebbsonlinecouple
!
!
!
!
  subroutine readsuewsout()
!
  use modulesuewsstebbscouple
!
  implicit none
!
  integer :: i,reason,icrop
!
!
!
  sout%timestep = 3600 ! 1hr, hard-coded for the test
!
!
!
  open(8,file='./SUEWS_output_res.csv',form='formatted')

  i = 0
  do while (.true.)
  read(8,*,iostat=reason)
  if( reason < 0 ) go to 333
  i = i + 1
  enddo
!
  333 continue
!
  sout%ntstep = i - 1
!
  allocate(sout%datetime(sout%ntstep))
  allocate(sout%hourmin(sout%ntstep))
  allocate(sout%Tair(   sout%ntstep))
  allocate(sout%Tsurf(  sout%ntstep))
  allocate(sout%Kwall(  sout%ntstep))
  allocate(sout%Kroof(  sout%ntstep))
  allocate(sout%ws(     sout%ntstep))
  allocate(sout%Lroof(  sout%ntstep))
  allocate(sout%Lwall(  sout%ntstep))
  allocate(sout%datetime_exch(sout%ntstep))
  allocate(sout%hourmin_exch(sout%ntstep))
  allocate(sout%Tair_exch(   sout%ntstep))
  allocate(sout%Tsurf_exch(  sout%ntstep))
  allocate(sout%Kwall_exch(  sout%ntstep))
  allocate(sout%Kroof_exch(  sout%ntstep))
  allocate(sout%ws_exch(     sout%ntstep))
  allocate(sout%Lroof_exch(  sout%ntstep))
  allocate(sout%Lwall_exch(  sout%ntstep))
!
!
!
  rewind(8)
!
  read(8,*)
!
  do i = 1,sout%ntstep,1
      read(8,*)sout%datetime(i),sout%hourmin(i),sout%Tair(i),sout%Tsurf(i), &
               sout%Kwall(i),sout%Kroof(i),sout%ws(i),sout%Lroof(i),sout%Lwall(i)
  enddo
!
  write(*,*)'    + SUEWS output profile'
  write(*,*)'    + Date            : ' , trim(sout%datetime(1)) , ' ' , trim(sout%hourmin(1)) , ' - ' , &
                                         trim(sout%datetime(sout%ntstep)), ' ' , trim(sout%hourmin(sout%ntstep))
  write(*,*)'    + Total data step : ' , sout%ntstep
!
  close(8)
!
!
!
!
  sout%ntskip = 12 
!
! This is the clock rate of original SUEWS by the STEBBS
! IF STEBBS uses 60 mins. output but SUEWS output 5 min. ntskip = 60/5 = 12 
!
!
  open(8,file='./SUEWS_output.csv',form='formatted')
!
  read(8,*)

  do i = 1,sout%ntstep*sout%ntskip,1
      if ( mod(i-1,sout%ntskip) == 0 ) then                          
          icrop = int((i-1)/sout%ntskip) + 1
          read(8,*)sout%datetime_exch(icrop),sout%hourmin_exch(icrop),  &
                   sout%Tair_exch(icrop),sout%Tsurf_exch(icrop),        &
                   sout%Kwall_exch(icrop),sout%Kroof_exch(icrop),       &
                   sout%ws_exch(icrop),sout%Lroof_exch(icrop),          &
                   sout%Lwall_exch(icrop)
      else
          read(8,*)
      endif
  enddo
!
  close(8)
!
!
!
  return
!
  end subroutine readsuewsout
!
!
!
!
!
  subroutine suewsstebbscouple(self)
!
  use modulestebbsprecision
  use modulestebbs, only : LBM, resolution
  use modulestebbsfunc, only : ext_conv_coeff
  use modulesuewsstebbscouple, only : sout,                                           &
                                      Tair_out, Tground_deep, Tsurf, density_air_out, &
                                      cp_air_out, Qsw_dn_extroof, Qsw_dn_extwall,     &
                                      Qlw_dn_extwall, Qlw_dn_extroof
!
  implicit none
!
  type(LBM) :: self
!
  integer :: tstep, i
  real(rprc)               :: Area, qfm_dom, qheat_dom, qcool_dom, qfb_hw_dom, qfb_dom_air, dom_temp, &
                              Qsw_transmitted_window, Qsw_absorbed_window, Qsw_absorbed_wallroof,     &
                              Qlw_net_extwallroof_to_outair, Qlw_net_extwindow_to_outair,             &
                              QStar, qinternal, qe_cool, qe_heat, QEC,                                &
                              Qconv_extwindow_to_outair, Qconv_extwallroof_to_outair, QH, QS,         &
                              Qcond_ground, Q_ventilation, QBAE, Q_waste, QWaste,                     &
                              temp, Textwallroof, Tintwallroof, Textwindow, Tintwindow, Tair_ind
!
  real(rprc),dimension(6)   :: bem_qf_1
  real(rprc),dimension(25)  :: energyEx
  character(len=256)              :: case
  character(len=256),dimension(4) :: fout
!
!
!
   case = self%case
   Area = self%Afootprint
!
!
!
  if ( self%flginit == 0 ) then
!
!     Output file
!
      fout(1) = 'Output_'//trim(case)//'.csv'; fout(2) = 'HeatFluxes_'//trim(case)//'.csv'
      fout(3) = 'EnergyBalance_'//trim(case)//'.csv'; fout(4) = 'Temp_'//trim(case)//'.csv'
!
      do i = 1, 4, 1
        open( i + 100*self%idLBM, file = trim(fout(i)), status = 'unknown' , form = 'formatted')
      enddo

      write(1+100*self%idLBM,*)',qheat_dom, qcool_dom, dom_tind, qfb_hw_dom, qfm_dom, qfb_dom_air'
      write(2+100*self%idLBM,*)',Qsw_transmitted_window, Qsw_absorbed_window, Qsw_absorbed_wallroof, '// &
                                'qmc, qir, qirw, qirf, qia, avr, qwc, qwic, qfc, qha, qwcon, qwicon, '// &
                   'qfcon, Qcond_ground, Qlw_net_extwallroof_to_outair, Qlw_net_extwindow_to_outair, '// &
                   'Qconv_extwallroof_to_outair, Qconv_extwindow_to_outair, qwaste, QS, QS_fabric, QS_air'
      write(3+100*self%idLBM,*)',QStar, QEC, QH, QS, QBAE, QWaste'
      write(4+100*self%idLBM,*)',Textwallroof, Tintwallroof, Textwindow, Tintwindow, Tair_ind'
!
  endif
!
!
!  Time integration start
!
  do tstep = 1, sout%ntstep, 1
!
       Tair_out        = sout%Tair(tstep) + 273.15
       Tground_deep    = 273.15 + 10.0
       Tsurf           = sout%Tsurf(tstep) + 273.15
       density_air_out = 1.225
       cp_air_out      = 1005.0
       Qsw_dn_extroof  = sout%Kroof(tstep)
       Qsw_dn_extwall  = sout%Kwall(tstep)
       Qlw_dn_extwall  = sout%Lwall(tstep)
       Qlw_dn_extroof  = sout%Lroof(tstep)
!
!
!
       self%h_o(1) = ext_conv_coeff(sout%ws_exch(tstep),sout%Tair_exch(tstep) - sout%Tsurf_exch(tstep))
       self%h_o(2) = ext_conv_coeff(sout%ws_exch(tstep),sout%Tair_exch(tstep) - sout%Tsurf_exch(tstep))
!
       call timeStepCalculation( self, Tair_out, Tground_deep, Tsurf,            &
                                 density_air_out, cp_air_out,                    &
                                 Qsw_dn_extroof, Qsw_dn_extwall,                 &
                                 Qlw_dn_extwall, Qlw_dn_extroof, sout%timestep,  &
                                 resolution                                      &
                               )
!
!
!       Handle return values
!
        bem_qf_1  = (/self%Qtotal_heating, self%Qtotal_cooling, self%EnergyExchanges(8),         &
                      self%Qtotal_water_tank, self%Qmetabolic_sensible, self%Qmetabolic_latent/)
        bem_qf_1  = bem_qf_1 / float(sout%timestep)
!
!       Metabolic sensible and latent heat
!
        qfm_dom   = bem_qf_1(5) + bem_qf_1(6)
!
!       Hourly heating load [W] and cooling load [W]
!
        qheat_dom = bem_qf_1(1)
        qcool_dom = bem_qf_1(2)
!
!       Hot water
!
        qfb_hw_dom = bem_qf_1(4)
!
!       Sensible heat to air [W]
!
        qfb_dom_air = 0 
!
        dom_temp = self%Tair_ind - 273.15 ! [K] to deg.
!
!
!
!        ?? = (/self%setTwater_tank, self%Twater_tank, self%Twater_vessel,                   &
!               self%Vwater_vessel, self%flowrate_water_supply, self%flowrate_water_drain/)
!
        energyEx = self%EnergyExchanges(:) / float(sout%timestep)
!
!       # calculate energy balance fluxes for a building, divided by footprint area [W m-2]
!       # 1) for net all wave radiation Q*
!       # knet from building = Qsw_transmitted_window + Qsw_absorbed_window + Qsw_absorbed_wallroof
!       # Lnet = -Qlw_net_extwallroof_to_outair - Qlw_net_extwindow_to_outair
!
        Qsw_transmitted_window = energyEx(1) / Area ! # transmitted solar radiation through windows [W m-2]
        Qsw_absorbed_window    = energyEx(2) / Area ! # absorbed solar radiation by windows [W m-2]
        Qsw_absorbed_wallroof  = energyEx(3) / Area ! #absorbed solar heat by walls [W m-2]
        Qlw_net_extwallroof_to_outair = energyEx(18) / Area ! # longwave radiation at external wall [W m-2]
        Qlw_net_extwindow_to_outair   = energyEx(19) / Area ! #longwave radiation at external windows [W m-2]
        QStar = Qsw_transmitted_window + Qsw_absorbed_window + Qsw_absorbed_wallroof  &
              - Qlw_net_extwallroof_to_outair - Qlw_net_extwindow_to_outair
!
!       # 2) sensible energy input into building (QEC)
!
        qinternal = (energyEx(8) + bem_qf_1(5))/Area ! #sensible internal appliance gain and sensible metabolism [W m-2]
        qe_cool = qcool_dom/self%coeff_performance_cooling/Area ! #energy use by cooling  [W m-2]
        qe_heat = qheat_dom/self%heating_efficiency_air/Area ! #energy use by heating [W m-2]
        QEC = qinternal + qe_cool + qe_heat ! # [W m-2] , Notice: energy use by hot water has not been added yet
!
!       # 3) sensible heat flux (QH)
!
        Qconv_extwindow_to_outair   = energyEx(21) / Area ! #convection at windows [W m-2]
        Qconv_extwallroof_to_outair = energyEx(20) / Area ! #convection at wall [W m-2]
        QH = Qconv_extwallroof_to_outair + Qconv_extwindow_to_outair ! #[W m-2]
!
!       # 4) storage heat flux (QS)
!
        qs = energyEx(23)/Area ! #heat storage/release by building fabric and indoor air  [W m-2]
        Qcond_ground = energyEx(17) / Area ! #conduction to external ground [W m-2], if assume ground floor is close to isolated, this flux should be close to 0
        QS = qs + Qcond_ground ! #[W m-2]
!
!       # 5) Building air exchange
!
        Q_ventilation = energyEx(9) / Area ! #ventilation and infiltration
        QBAE = -Q_ventilation ! #[W m-2]
!
!       # 6) Waste heat by mechani cooling (waste heat by heating is released into indoor, so not in this part)
!
        q_waste = energyEx(22)/Area ! # waste heat from cooling  include energy consumption
        QWaste = q_waste    ! #[W m-2]
!
!       Temperature output
!
        Textwallroof = self%Textwallroof ! # external surface temperature of wall [K]
        Tintwallroof = self%Tintwallroof ! # internal surface temperature of wall [K]
        Textwindow   = self%Textwindow   ! # external surface temperature of window [K]
        Tintwindow   = self%Tintwindow   ! # internal surface temperature of window [K]
        Tair_ind     = self%Tair_ind     ! # Indoor air temperature [K]
!
!
!
        write(1+100*self%idLBM,'(a19,1x,6(",",f10.5))' ) trim(sout%datetime(tstep))//' '//trim(sout%hourmin(tstep)),  &
                                                         qheat_dom, qcool_dom, dom_temp, qfb_hw_dom, qfm_dom, qfb_dom_air
        write(2+100*self%idLBM,'(a19,1x,25(",",f15.5))') trim(sout%datetime(tstep))//' '//trim(sout%hourmin(tstep)), &
                                      ( energyEx(i)/Area, i = 1, 25, 1)
        write(3+100*self%idLBM,'(a19,1x,6(",",f15.5))' ) trim(sout%datetime(tstep))//' '//trim(sout%hourmin(tstep)), QStar, QEC,  &
                                                         QH, QS, QBAE, QWaste
        write(4+100*self%idLBM,'(a19,1x,5(",",f15.5))' ) trim(sout%datetime(tstep))//' '//trim(sout%hourmin(tstep)), Textwallroof, &
                                                         Tintwallroof, Textwindow, Tintwindow, Tair_ind
!
!
!
  enddo
!
!
! Store the STEBBS output to building type variables
!
  self%qfm_dom     = qfm_dom
  self%qheat_dom   = qheat_dom
  self%qcool_dom   = qcool_dom
  self%qfb_hw_dom  = qfb_hw_dom
  self%dom_temp    = dom_temp
  self%QStar       = QStar
  self%QEC         = QEC
  self%QH          = QH
  self%QS          = QS
  self%QBAE        = QBAE
  self%QWaste      = QWaste
!
!
!
  self%flginit = 1
!
  return
!
  end subroutine suewsstebbscouple
!
!
!
!
  subroutine timeStepCalculation( self, Tair_out, Tground_deep, Tsurf, &
                                  density_air_out, cp_air_out,         &
                                  Qsw_dn_extroof, Qsw_dn_extwall,      &
                                  Qlw_dn_extwall, Qlw_dn_extroof,      &
                                  timestep, resolution                 &
                                 )
!
  use modulestebbsprecision
  use modulestebbs, only: LBM
!
  implicit none
!
  integer          :: timestep, resolution
  real(rprc)       :: Tair_out, Tground_deep, Tsurf, density_air_out, &
                      cp_air_out, Qsw_dn_extroof, Qsw_dn_extwall,     &
                      Qlw_dn_extwall, Qlw_dn_extroof
!
  type(LBM) :: self
!
!
! ReinitialiseHC
!
  self%Qtotal_heating = 0.0
  self%Qtotal_cooling = 0.0
  self%Qtotal_water_tank = 0.0
  self%qhwtDrain = 0.0
!
!
!
  call tstep(                                                                           &
  Tair_out, Tground_deep, Tsurf,                                                        &
  density_air_out, cp_air_out,                                                          &
  Qsw_dn_extroof, Qsw_dn_extwall,                                                       &
  Qlw_dn_extwall, Qlw_dn_extroof,                                                       &
  self%wiTAR(1), self%wiTAR(2), self%wiTAR(3),                                          &
  self%waTAR(1), self%waTAR(2), self%waTAR(3),                                          &
  self%Qtotal_heating, self%Qtotal_cooling,                                             &
  self%height_building, self%ratio_window_wall, self%thickness_wallroof,                &
  self%thickness_groundfloor, self%depth_ground, self%thickness_window,                 &
  self%conv_coeff_intwallroof, self%conv_coeff_indoormass,                              &
  self%conv_coeff_intgroundfloor, self%conv_coeff_intwindow,                            &
!  self%conv_coeff_extwallroof, self%conv_coeff_extwindow,                               &
  self%h_o(1), self%h_o(2),                                                             &
  self%conductivity_wallroof, self%conductivity_groundfloor,                            &
  self%conductivity_window, self%conductivity_ground,                                   &
  self%density_wallroof, self%density_groundfloor,                                      &
  self%density_window, self%density_indoormass, self%density_air_ind,                   &
  self%cp_wallroof, self%cp_groundfloor, self%cp_window,                                &
  self%cp_indoormass, self%cp_air_ind,                                                  &
  self%emissivity_extwallroof, self%emissivity_intwallroof, self%emissivity_indoormass, &
  self%emissivity_extwindow, self%emissivity_intwindow,                                 &
  self%windowTransmissivity, self%windowAbsorbtivity, self%windowReflectivity,          &
  self%wallTransmisivity, self%wallAbsorbtivity, self%wallReflectivity,                 &
  self%BVF_extwall, self%GVF_extwall, self%SVF_extwall,                                 &
  self%occupants, self%metabolic_rate, self%ratio_metabolic_latent_sensible,            &
  self%appliance_power_rating, self%appliance_usage_factor,                             &
  self%maxheatingpower_air, self%heating_efficiency_air,                                &
  self%maxcoolingpower_air, self%coeff_performance_cooling,                             &
  self%Vair_ind, self%ventilation_rate, self%Awallroof,                                 &
  self%Vwallroof, self%Afootprint, self%Vgroundfloor,                                   &
  self%Awindow, self%Vwindow, self%Vindoormass, self%Aindoormass,                       &
  self%Tair_ind , self%Tindoormass, self%Tintwallroof, self%Textwallroof,               &
  self%Tintwindow, self%Textwindow, self%Tintgroundfloor, self%Textgroundfloor,         &
  self%Ts,                                                                              &
!  self%Ts(1), self%Ts(2),                                                               &
  self%appliance_totalnumber, timestep, resolution,                                     &
  self%Qtotal_water_tank, self%Twater_tank, self%Tintwall_tank,                         &
  self%Textwall_tank, self%thickness_tankwall, self%Tincomingwater_tank,                &
  self%Vwater_tank, self%Asurf_tank, self%Vwall_tank, self%setTwater_tank,              &
  self%Twater_vessel, self%Tintwall_vessel, self%Textwall_vessel,                       &
  self%thickness_wall_vessel, self%Vwater_vessel, self%Awater_vessel,                   &
  self%Vwall_vessel, self%flowrate_water_supply, self%flowrate_water_drain,             &
  self%cp_water, self%cp_wall_tank, self%cp_wall_vessel,                                &
  self%density_water, self%density_wall_tank, self%density_wall_vessel,                 &
  self%BVF_tank, self%MVF_tank, self%conductivity_wall_tank,                            &
  self%conv_coeff_intwall_tank, self%conv_coeff_extwall_tank,                           &
  self%emissivity_extwall_tank, self%conductivity_wall_vessel,                          &
  self%conv_coeff_intwall_vessel, self%conv_coeff_extwall_vessel,                       &
  self%emissivity_extwall_vessel, self%HeatingPower_DHW,                                &
  self%heating_efficiency_water, self%minVwater_vessel,                                 &
  self%weighting_factor_heatcapacity_wallroof,                                          &
!
! Output only variables
!
  self%EnergyExchanges(1),  & !Qsw_transmitted_window_tstepTotal
  self%EnergyExchanges(2),  & !Qsw_absorbed_window_tstepTotal
  self%EnergyExchanges(3),  & !Qsw_absorbed_wallroof_tstepTotal
  self%EnergyExchanges(4),  & !Qconv_indair_to_indoormass_tstepTotal
  self%EnergyExchanges(5),  & !Qlw_net_intwallroof_to_allotherindoorsurfaces_tstepTotal
  self%EnergyExchanges(6),  & !Qlw_net_intwindow_to_allotherindoorsurfaces_tstepTotal
  self%EnergyExchanges(7),  & !Qlw_net_intgroundfloor_to_allotherindoorsurfaces_tstepTotal
  self%EnergyExchanges(8),  & !Q_appliance_tstepTotal
  self%EnergyExchanges(9),  & !Q_ventilation_tstepTotal
  self%EnergyExchanges(10), & !Qconv_indair_to_intwallroof_tstepTotal
  self%EnergyExchanges(11), & !Qconv_indair_to_intwindow_tstepTotal
  self%EnergyExchanges(12), & !Qconv_indair_to_intgroundfloor_tstepTotal
  self%EnergyExchanges(13), & !Qloss_efficiency_heating_air_tstepTotal
  self%EnergyExchanges(14), & !Qcond_wallroof_tstepTotal
  self%EnergyExchanges(15), & !Qcond_window_tstepTotal
  self%EnergyExchanges(16), & !Qcond_groundfloor_tstepTotal
  self%EnergyExchanges(17), & !Qcond_ground_tstepTotal
  self%EnergyExchanges(18), & !Qlw_net_extwallroof_to_outair_tstepTotal
  self%EnergyExchanges(19), & !Qlw_net_extwindow_to_outair_tstepTotal
  self%EnergyExchanges(20), & !Qconv_extwallroof_to_outair_tstepTotal
  self%EnergyExchanges(21), & !Qconv_extwindow_to_outair_tstepTotal
  self%EnergyExchanges(22), & !q_cooling_timestepTotal
  self%EnergyExchanges(23), & !QS_tstepTotal
  self%EnergyExchanges(24), & !QS_fabric_tstepTotal
  self%EnergyExchanges(25), & !QS_air_tstepTotal
  self%qhwtDrain,           & !Qloss_drain
  self%Qmetabolic_sensible, & !qsensible_timestepTotal
  self%Qmetabolic_latent)     !qlatent_timestepTotal
!
!
!
  return
!
  end subroutine timeStepCalculation
!
!
!
!
!
  subroutine tstep(                                                      &
  Tair_out, Tground_deep, Tsurf,                                         &
  density_air_out, cp_air_out,                                           &
  Qsw_dn_extroof, Qsw_dn_extwall,                                        &
  Qlw_dn_extwall, Qlw_dn_extroof,                                        &
  winT, winA, winR, walT, walA, walR,                                    &
  Qtotal_heating, Qtotal_cooling,                                        & !IO
  height_building, ratio_window_wall, thickness_wallroof,                &
  thickness_groundfloor, depth_ground, thickness_window,                 &
  conv_coeff_intwallroof, conv_coeff_indoormass,                         &
  conv_coeff_intgroundfloor, conv_coeff_intwindow,                       &
  conv_coeff_extwallroof, conv_coeff_extwindow,                          &
  conductivity_wallroof, conductivity_groundfloor,                       &
  conductivity_window, conductivity_ground,                              &
  density_wallroof, density_groundfloor,                                 &
  density_window, density_indoormass, density_air_ind,                   &
  cp_wallroof, cp_groundfloor, cp_window, cp_indoormass, cp_air_ind,     &
  emissivity_extwallroof, emissivity_intwallroof, emissivity_indoormass, &
  emissivity_extwindow, emissivity_intwindow,                            &
  windowTransmissivity, windowAbsorbtivity, windowReflectivity,          &
  wallTransmisivity, wallAbsorbtivity, wallReflectivity,                 &
  BVF_extwall, GVF_extwall, SVF_extwall,                                 &
  occupants, metabolic_rate, ratio_metabolic_latent_sensible,            &
  appliance_power_rating, appliance_usage_factor,                        &
  maxheatingpower_air, heating_efficiency_air,                           &
  maxcoolingpower_air, coeff_performance_cooling,                        &
  Vair_ind, ventilation_rate, Awallroof,                                 &
  Vwallroof, Afootprint, Vgroundfloor,                                   &
  Awindow, Vwindow, Vindoormass, Aindoormass,                            &
  Tair_ind , Tindoormass, Tintwallroof, Textwallroof,                    & !IO
  Tintwindow, Textwindow, Tintgroundfloor, Textgroundfloor,              & !IO
  Ts,                                                                    & !IO
!  Ts(1), Ts(2),                                                          &
  appliance_totalnumber, timestep, resolution,                           &
  Qtotal_water_tank, Twater_tank, Tintwall_tank,                         & !IO
  Textwall_tank, thickness_tankwall, Tincomingwater_tank,                & !IO
  Vwater_tank, Asurf_tank, Vwall_tank, setTwater_tank,                   & !IO
  Twater_vessel, Tintwall_vessel, Textwall_vessel,                       & !IO
  thickness_wall_vessel, Vwater_vessel, Awater_vessel,                   & !IO
  Vwall_vessel, flowrate_water_supply, flowrate_water_drain,             &
  cp_water, cp_wall_tank, cp_wall_vessel,                                &
  density_water, density_wall_tank, density_wall_vessel,                 &
  BVF_tank, MVF_tank, conductivity_wall_tank,                            &
  conv_coeff_intwall_tank, conv_coeff_extwall_tank,                      &
  emissivity_extwall_tank, conductivity_wall_vessel,                     &
  conv_coeff_intwall_vessel, conv_coeff_extwall_vessel,                  &
  emissivity_extwall_vessel, maxheatingpower_water,                      &
  heating_efficiency_water, minVwater_vessel,                            &
  weighting_factor_heatcapacity_wallroof,                                &
!
! Output only variables
!
  Qsw_transmitted_window_tstepTotal,                                     & !EE(1)
  Qsw_absorbed_window_tstepTotal,                                        & !EE(2)
  Qsw_absorbed_wallroof_tstepTotal,                                      & !EE(3)
  Qconv_indair_to_indoormass_tstepTotal,                                 & !EE(4)
  Qlw_net_intwallroof_to_allotherindoorsurfaces_tstepTotal,              & !EE(5)
  Qlw_net_intwindow_to_allotherindoorsurfaces_tstepTotal,                & !EE(6)
  Qlw_net_intgroundfloor_to_allotherindoorsurfaces_tstepTotal,           & !EE(7)
  Q_appliance_tstepTotal,                                                & !EE(8)
  Q_ventilation_tstepTotal, Qconv_indair_to_intwallroof_tstepTotal,      & !EE(9),(10)
  Qconv_indair_to_intwindow_tstepTotal,                                  & !EE(11)
  Qconv_indair_to_intgroundfloor_tstepTotal,                             & !EE(12)
  Qloss_efficiency_heating_air_tstepTotal,                               & !EE(13)
  Qcond_wallroof_tstepTotal, Qcond_window_tstepTotal,                    & !EE(14),(15)
  Qcond_groundfloor_tstepTotal, Qcond_ground_tstepTotal,                 & !EE(16),(17)
  Qlw_net_extwallroof_to_outair_tstepTotal,                              & !EE(18)
  Qlw_net_extwindow_to_outair_tstepTotal,                                & !EE(19)
  Qconv_extwallroof_to_outair_tstepTotal,                                & !EE(20)
  Qconv_extwindow_to_outair_tstepTotal,                                  & !EE(21)
  q_cooling_timestepTotal,                                               & !EE(22)
  QS_tstepTotal, QS_fabric_tstepTotal, QS_air_tstepTotal,                & !EE(23,24,25)
  Qloss_drain,                                                           & !qhwtDrain
  qsensible_timestepTotal, qlatent_timestepTotal)                          !Qmetabolic_sensible, Qmetabolic_latent
!
  use modulestebbsprecision
  use modulestebbsfunc
!
  implicit none
!
  integer          :: i
  real(rprc)       :: Tair_out, Tground_deep, Tsurf,                  &
                      density_air_out, cp_air_out, Qsw_dn_extroof,    &
                      Qsw_dn_extwall, Qlw_dn_extwall, Qlw_dn_extroof
!    /*** DOMESTIC HOT WATER ***/
  real(rprc)       :: Twater_tank,   & ! // Water temperature (K) in Hot Water Tank
                      Tintwall_tank, & ! // Hot water tank internal wall temperature (K)
                      Textwall_tank    ! //Hot water tank external wall temperature (K)
  real(rprc)       :: dTwater_tank = 0.0, dTintwall_tank = 0.0, dTextwall_tank = 0.0
  real(rprc)       :: thickness_tankwall ! //Hot water tank wall thickness
!
  real(rprc)       :: Tincomingwater_tank ! // Water temperature (K) of Water coming into the Water Tank
  real(rprc)       :: Vwater_tank, & ! // Volume of Water in Hot Water Tank (m^3)
                      Asurf_tank,  & ! // Surface Area of Hot Water Tank(m^2)
                      Vwall_tank,  & ! // Wall volume of Hot Water Tank(m^2)
                      setTwater_tank ! // Water Tank setpoint temperature (K)
!
  real(rprc)       :: Twater_vessel,   & ! // Water temperature (K) of water held in use in Building
                      Tintwall_vessel, & ! // Hot water tank internal wall temperature (K)
                      Textwall_vessel    ! // Hot water tank external wall temperature (K)
  real(rprc)       :: dTwater_vessel = 0.0, dTintwall_vessel = 0.0, dTextwall_vessel = 0.0
  real(rprc)       :: thickness_wall_vessel ! //DHW vessels wall thickness
!
  real(rprc)       :: Vwater_vessel            ! // Volume of water held in use in building
  real(rprc)       :: dVwater_vessel = 0.0     ! // Change in volume of Domestic Hot Water held in use in building
  real(rprc)       :: Awater_vessel,         & ! // Surface Area of Hot Water in Vessels in Building
                      Vwall_vessel,          & ! // Wall volume of Hot water Vessels in Building
                      flowrate_water_supply, & ! // Hot Water Flow Rate in m^3 / s
                      flowrate_water_drain     ! // Draining of Domestic Hot Water held in building
!
  real(rprc)       :: cp_water,      & ! //Specific Heat Capacity of Domestic Hot Water
                      cp_wall_tank,  & ! //Specific Heat Capacity of Hot Water Tank wall
                      cp_wall_vessel   ! //Specific Heat Capacity of Vessels containing DHW in use in Building
!
  real(rprc)       :: density_water,     & ! //Density of water
                      density_wall_tank, & ! //Density of hot water tank wall
                      density_wall_vessel  ! //Density of vessels containing DHW in use in buildings
!
  real(rprc)       :: BVF_tank, & ! //water tank - building wall view factor
                      MVF_tank    ! //water tank - building internal mass view factor
!
  real(rprc)       :: conductivity_wall_tank,    & ! // Effective Wall conductivity of the Hot Water Tank
                      conv_coeff_intwall_tank,   & ! // Effective Internal Wall convection coefficient of the Hot Water Tank
                      conv_coeff_extwall_tank,   & ! // Effective External Wall convection coefficient of the Hot Water Tank
                      emissivity_extwall_tank,   & ! // Effective External Wall emissivity of the Hot Water Tank
                      conductivity_wall_vessel,  & ! // Effective Conductivity of vessels containing DHW in use in Building.
                      conv_coeff_intwall_vessel, & ! // Effective Internal Wall convection coefficient of the Vessels holding DHW in use in Building
                      conv_coeff_extwall_vessel, & ! // Effective Enternal Wall convection coefficient of the Vessels holding DHW in use in Building
                      emissivity_extwall_vessel    ! // Effective External Wall emissivity of hot water being used within building
!    /** NOTE THAT LATENT HEAT FLUX RELATING TO HOT WATER USE CURRENTLY NOT IMPLEMENTED **/
!
  real(rprc)       :: maxheatingpower_water,  &
                      heating_efficiency_water
!    /*** END DOMESTIC HOT WATER ***/
!
  real(rprc)       :: winT, & ! // window transmisivity
                      winA, & ! // window absorptivity
                      winR, & ! // window reflectivity
                      walT, & ! // wall transmisivity
                      walA, & ! // wall absorptivity
                      walR    ! // wall reflectivity
!
  real(rprc)       :: Qtotal_heating, & ! // currently only sensible but this needs to be  split into sensible and latent heat components
                      Qtotal_cooling    ! // currently only sensible but this needs to be  split into sensible and latent heat components
!
  real(rprc)       :: height_building, ratio_window_wall, &
                      thickness_wallroof, thickness_groundfloor, depth_ground, thickness_window, & 
!    //float height_building, width, depth, ratio_window_wall, thickness_wallroof, thickness_groundfloor, depth_ground, thickness_window;
                      conv_coeff_intwallroof, conv_coeff_indoormass,         &
                      conv_coeff_intgroundfloor, conv_coeff_intwindow,       &
                      conv_coeff_extwallroof, conv_coeff_extwindow,          &
                      conductivity_wallroof, conductivity_groundfloor,       &
                      conductivity_window, conductivity_ground,              &
                      density_wallroof, density_groundfloor, density_window, &
                      density_indoormass, density_air_ind,                   &
                      cp_wallroof, cp_groundfloor, cp_window,                &
                      cp_indoormass, cp_air_ind,                             &
                      emissivity_extwallroof, emissivity_intwallroof,        &
                      emissivity_indoormass, emissivity_extwindow, emissivity_intwindow, &
                      windowTransmissivity, windowAbsorbtivity, windowReflectivity,      &
                      wallTransmisivity, wallAbsorbtivity, wallReflectivity,             &
                      BVF_extwall, GVF_extwall, SVF_extwall
!
  real(rprc)       :: occupants
!
  real(rprc)       :: metabolic_rate, ratio_metabolic_latent_sensible, &
                      appliance_power_rating
  integer          :: appliance_totalnumber
  real(rprc)       :: appliance_usage_factor,                         &
                      maxheatingpower_air, heating_efficiency_air,    &
                      maxcoolingpower_air, coeff_performance_cooling, &
                      Vair_ind, ventilation_rate, & ! //Fixed at begining to have no natural ventilation. Given in units of volume of air per hour
                      Awallroof, Vwallroof,                           &
                      Afootprint, Vgroundfloor,                       &
                      Awindow, Vwindow,                               &
                      Vindoormass, Aindoormass ! //Assumed internal mass as a cube
!
  real(rprc)       :: Tair_ind, Tindoormass, Tintwallroof, Textwallroof,         &
                      Tintwindow, Textwindow, Tintgroundfloor, Textgroundfloor
  real(rprc)       :: dTair_ind = 0.0, dTindoormass = 0.0, dTintwallroof = 0.0,  &
                      dTextwallroof = 0.0, dTintwindow = 0.0, dTextwindow = 0.0, &
                      dTintgroundfloor = 0.0, dTextgroundfloor = 0.0
!
!    /*** DOMESTIC HOT WATER ***/
  real(rprc)       :: Qconv_water_to_inttankwall = 0.0,          & ! // heat flux to internal wall of hot water tank
                      Qconv_exttankwall_to_indair = 0.0,         & ! // convective heat flux to external wall of hot water tank
                      Qlw_net_exttankwall_to_intwallroof = 0.0,  & ! 
                      Qlw_net_exttankwall_to_indoormass = 0.0,   & ! // radiative heat flux to external wall of hot water tank
                      Qcond_tankwall = 0.0,                      & ! // heat flux through wall of hot water tank
                      Qtotal_water_tank,                         & ! // total heat input into water of hot water tank over simulation, hence do not equate to zero
                      Qconv_water_to_intvesselwall = 0.0,        & ! // heat flux to internal wall of vessels holding DHW in use in building
                      Qcond_vesselwall = 0.0,                    & ! // heat flux through wall of vessels holding DHW in use in building
                      Qconv_extvesselwall_to_indair = 0.0,       & ! // convective heat flux to external wall of vessels holding DHW in use in building
                      Qlw_net_extvesselwall_to_wallroof = 0.0,   & 
                      Qlw_net_extvesselwall_to_indoormass = 0.0, & ! // radiative heat flux to external wall of vessels holding DHW in use in building
!                      Qloss_drain = 0.0,                         & ! //Heat loss as water held in use in building drains to sewer
                      Qloss_efficiency_heating_water = 0.0 ! // additional heat release from efficieny losses/gains of heating hot water

  real(rprc), intent(out) :: Qloss_drain           ! // Heat loss as water held in use in building drains to sewer
  real(rprc)       :: Qtotal_net_water_tank = 0.0, Qtotal_net_intwall_tank = 0.0,      &
                      Qtotal_net_extwall_tank = 0.0, Qtotal_net_water_vessel = 0.0,    &
                      Qtotal_net_intwall_vessel = 0.0, Qtotal_net_extwall_vessel = 0.0
  real(rprc)       :: qhwt_timestep=0.0
!
  real(rprc)       :: VARatio_water_vessel = 0.0
  real(rprc)       :: minVwater_vessel
  real(rprc)       :: weighting_factor_heatcapacity_wallroof
!    /************************************************************/
!    /*** END DOMESTIC HOT WATER ***/
!
  real(rprc),dimension(2) :: Ts ! //Heating and Cooling setpoint temperature (K)s, respectively
  real(rprc),dimension(2) :: Qm ! //Metabolic heat, sensible(1) and latent(2)
!
  integer :: timestep, resolution
!
  real(rprc)       :: Qf_ground_timestep = 0.0, &
                      q_heating_timestep = 0.0, &
                      q_cooling_timestep = 0.0
!
  real(rprc)       :: Qsw_transmitted_window = 0.0, Qsw_absorbed_window = 0.0, Qsw_absorbed_wallroof = 0.0,  &
                      Qconv_indair_to_indoormass = 0.0, Qlw_net_intwallroof_to_allotherindoorsurfaces = 0.0, &
                      Qlw_net_intwindow_to_allotherindoorsurfaces = 0.0, Qlw_net_intgroundfloor_to_allotherindoorsurfaces = 0.0
  real(rprc)       :: Q_appliance = 0.0, Q_ventilation = 0.0, Qconv_indair_to_intwallroof = 0.0,  &
                      Qconv_indair_to_intwindow = 0.0, Qconv_indair_to_intgroundfloor = 0.0
  real(rprc)       :: Qloss_efficiency_heating_air = 0.0, Qcond_wallroof = 0.0, Qcond_window = 0.0, &
                      Qcond_groundfloor = 0.0, Qcond_ground = 0.0
  real(rprc)       :: Qlw_net_extwallroof_to_outair = 0.0, Qlw_net_extwindow_to_outair = 0.0, &
                      Qconv_extwallroof_to_outair = 0.0, Qconv_extwindow_to_outair = 0.0
  real(rprc)       :: QS_total = 0.0, QS_fabric = 0.0, QS_air = 0.0
!    //STS - 31/07/2018: Added parameters to return surface fluxes for each timestep.
!    //This allows model user to calculate QfB and Qs externally to the model
!    //and provide a more flexible interface to surface energy balance models.
!
! Output vars
  real(rprc),intent(inout)       :: Qsw_transmitted_window_tstepTotal,                        &
                                    Qsw_absorbed_window_tstepTotal,                           &
                                    Qsw_absorbed_wallroof_tstepTotal,                         &
                                    Qconv_indair_to_indoormass_tstepTotal,                    &
                                    Qlw_net_intwallroof_to_allotherindoorsurfaces_tstepTotal, &
                                    Qlw_net_intwindow_to_allotherindoorsurfaces_tstepTotal,   & 
                                    Qlw_net_intgroundfloor_to_allotherindoorsurfaces_tstepTotal
  real(rprc),intent(inout)       :: Q_appliance_tstepTotal,                   &
                                    Q_ventilation_tstepTotal,                 &
                                    Qconv_indair_to_intwallroof_tstepTotal,   &
                                    Qconv_indair_to_intwindow_tstepTotal,     &
                                    Qconv_indair_to_intgroundfloor_tstepTotal
  real(rprc),intent(inout)       :: Qloss_efficiency_heating_air_tstepTotal,  &
                                    Qcond_wallroof_tstepTotal,                &
                                    Qcond_window_tstepTotal,                  &
                                    Qcond_groundfloor_tstepTotal,             &
                                    Qcond_ground_tstepTotal
  real(rprc),intent(inout)       :: Qlw_net_extwallroof_to_outair_tstepTotal, &
                                    Qlw_net_extwindow_to_outair_tstepTotal,   &
                                    Qconv_extwallroof_to_outair_tstepTotal,   &
                                    Qconv_extwindow_to_outair_tstepTotal
  real(rprc),intent(inout)       :: q_cooling_timestepTotal, &
                                    qsensible_timestepTotal, &
                                    qlatent_timestepTotal
  real(rprc),intent(inout)       :: QS_tstepTotal, QS_fabric_tstepTotal, QS_air_tstepTotal
  real(rprc) :: Qmetabolic_sensible = 0.0, Qmetabolic_latent = 0.0
!
  real(rprc)       :: Qtotal_net_indoormass = 0.0, Qtotal_net_indair = 0.0,        &
                      Qtotal_net_intwallroof = 0.0, Qtotal_net_extwallroof = 0.0,  &
                      Qtotal_net_intwindow = 0.0, Qtotal_net_extwindow = 0.0,      &
                      Qtotal_net_intgroundfloor = 0.0, Qtotal_net_extgroundfloor = 0.0
!
!
!
! Output cleaning
!
  Qsw_transmitted_window_tstepTotal = 0.0;                        &
  Qsw_absorbed_window_tstepTotal = 0.0;                           &
  Qsw_absorbed_wallroof_tstepTotal = 0.0;                         &
  Qconv_indair_to_indoormass_tstepTotal = 0.0;                    &
  Qlw_net_intwallroof_to_allotherindoorsurfaces_tstepTotal = 0.0; &
  Qlw_net_intwindow_to_allotherindoorsurfaces_tstepTotal = 0.0;   &
  Qlw_net_intgroundfloor_to_allotherindoorsurfaces_tstepTotal = 0.0
  Q_appliance_tstepTotal = 0.0;                   &
  Q_ventilation_tstepTotal = 0.0;                 &
  Qconv_indair_to_intwallroof_tstepTotal = 0.0;   &
  Qconv_indair_to_intwindow_tstepTotal = 0.0;     &
  Qconv_indair_to_intgroundfloor_tstepTotal = 0.0
  Qloss_efficiency_heating_air_tstepTotal = 0.0;  &
  Qcond_wallroof_tstepTotal = 0.0;                &
  Qcond_window_tstepTotal = 0.0;                  &
  Qcond_groundfloor_tstepTotal = 0.0;             &
  Qcond_ground_tstepTotal = 0.0
  Qlw_net_extwallroof_to_outair_tstepTotal = 0.0; &
  Qlw_net_extwindow_to_outair_tstepTotal = 0.0;   &
  Qconv_extwallroof_to_outair_tstepTotal = 0.0;   &
  Qconv_extwindow_to_outair_tstepTotal = 0.0
  q_cooling_timestepTotal = 0.0; &
  qsensible_timestepTotal = 0.0; &
  qlatent_timestepTotal = 0.0
  QS_tstepTotal = 0.0; QS_fabric_tstepTotal = 0.0; QS_air_tstepTotal = 0.0
!
!
!
!
! Simulation starts
!
!    //Used to recalculate Area of DHW in use
    if (Awater_vessel > 0.0) then
        VARatio_water_vessel = Vwater_vessel / Awater_vessel
    endif
!
!    if ((timestep % resolution) == 0) {
!        for (int i=0;i<(timestep/resolution);i++) {
    if ( mod(timestep, resolution) == 0 ) then
        looptime: do i = 1, int(timestep/resolution), 1
!
            Qsw_transmitted_window = windowInsolation(Qsw_dn_extwall, winT, Awindow)
            Qsw_absorbed_window = windowInsolation(Qsw_dn_extwall, winA, Awindow)
!            // Awallroof excludes windows and includes floor area
            Qsw_absorbed_wallroof =                                         &
            wallInsolation(Qsw_dn_extwall, walA, Awallroof - Afootprint) +  &
            wallInsolation(Qsw_dn_extroof, walA, Afootprint) !//separate the wall and roof
!//            printf("Qconv_indair_to_indoormass: %f  Tindoormass: %f  Tair_ind: %f  ", Qconv_indair_to_indoormass, Tindoormass, Tair_ind);
            Qconv_indair_to_indoormass = internalConvectionHeatTransfer(conv_coeff_indoormass, Aindoormass, Tindoormass, Tair_ind)
!//            printf("Qconv_indair_to_indoormass: %f  Tindoormass: %f  Tair_ind: %f \n", Qconv_indair_to_indoormass, Tindoormass, Tair_ind);
            Qlw_net_intwallroof_to_allotherindoorsurfaces = indoorRadiativeHeatTransfer() ! //  for wall internal radiative exchange
            Qlw_net_intwindow_to_allotherindoorsurfaces = Qlw_net_intwallroof_to_allotherindoorsurfaces ! //  for window internal radiative exchange - TODO: currently no distinction in internal radiative exchanges
            Qlw_net_intgroundfloor_to_allotherindoorsurfaces = Qlw_net_intwallroof_to_allotherindoorsurfaces ! //  for ground floor internal radiative exchange - TODO: currently no distinction in internal radiative exchanges
!//
            Q_appliance =                                                                                  &
            internalApplianceGains(appliance_power_rating, appliance_usage_factor, appliance_totalnumber)
            Q_ventilation =                                                                                &
            ventilationHeatTransfer(density_air_ind, cp_air_ind, ventilation_rate, Tair_out, Tair_ind)
            Qconv_indair_to_intwallroof =                                                                  &
            indoorConvectionHeatTransfer(conv_coeff_intwallroof, Awallroof, Tintwallroof, Tair_ind)
            Qconv_indair_to_intwindow =                                                                    &
            indoorConvectionHeatTransfer(conv_coeff_intwindow, Awindow, Tintwindow, Tair_ind)
            Qconv_indair_to_intgroundfloor =                                                               &
            indoorConvectionHeatTransfer(conv_coeff_intgroundfloor, Afootprint, Tintgroundfloor, Tair_ind)
!//
            q_heating_timestep = heating(Ts(1), Tair_ind, heating_efficiency_air, maxheatingpower_air)
            q_cooling_timestep = cooling(Ts(2), Tair_ind, coeff_performance_cooling, maxcoolingpower_air)

            !internalOccupancyGains(occupants, metabolic_rate, ratio_metabolic_latent_sensible, Qmetabolic_sensible, Qmetabolic_latent)
            Qm = internalOccupancyGains(occupants, metabolic_rate, ratio_metabolic_latent_sensible)
!
            Qmetabolic_sensible = Qm(1) 
            Qmetabolic_latent   = Qm(2)
!
            Qloss_efficiency_heating_air =                                                                 &
            additionalSystemHeatingEnergy(q_heating_timestep, heating_efficiency_air)
            Qcond_wallroof =                                                                               &
            wallConduction(conductivity_wallroof, Awallroof, Tintwallroof, Textwallroof, thickness_wallroof)
            Qcond_window =                                                                                 &
            windowConduction(conductivity_window, Awindow, Tintwindow, Textwindow, thickness_window)
            Qcond_groundfloor =                                                                            &
            wallConduction(conductivity_groundfloor, Afootprint, Tintgroundfloor, Textgroundfloor, thickness_groundfloor)
            Qcond_ground =                                                                                 &
            wallConduction(conductivity_ground, Afootprint, Textgroundfloor, Tground_deep, depth_ground) ! ; // conduction from ground floor to external ground - depth of the ground can be set (depth_ground) and ground temperature should be determined accordingly..
            ! //add the longwave radiation between sky and external envelope; Assume surrounding surface and ground surface temperature is the same (e.g. assumed outdoor air temperature as in EnergyPlus)
            ! // Qlw_net_extwallroof_to_outair = outdoorRadiativeHeatTransfer(BVF_extwall, Awallroof, emissivity_extwallroof, Textwallroof, Tsurf)+outdoorRadiativeHeatTransfer(SVF_extwall, Awallroof, emissivity_extwallroof, Textwallroof, Tsky);
            ! // Qlw_net_extwindow_to_outair = outdoorRadiativeHeatTransfer(BVF_extwall, Awindow, emissivity_extwindow, Textwindow, Tsurf)+outdoorRadiativeHeatTransfer(SVF_extwall, Awindow, emissivity_extwindow, Textwindow, Tsky);
            ! // call outdoorRadiativeHeatTransfer with LW instead of temp
            Qlw_net_extwallroof_to_outair =                                                                &
            lwoutdoorRadiativeHeatTransfer                                                                 &
            (Awallroof, emissivity_extwallroof, Textwallroof,                                              &
             ((Qlw_dn_extwall * (Awallroof - Afootprint)) + (Qlw_dn_extroof * Afootprint)) / Awallroof)
            Qlw_net_extwindow_to_outair = lwoutdoorRadiativeHeatTransfer(Awindow, emissivity_extwindow, Textwindow, Qlw_dn_extwall)
            Qconv_extwallroof_to_outair = outdoorConvectionHeatTransfer(conv_coeff_extwallroof, Awallroof, Textwallroof, Tair_out)
            Qconv_extwindow_to_outair = outdoorConvectionHeatTransfer(conv_coeff_extwindow, Awindow, Textwindow, Tair_out)

!
!            /**************************/
!            /*** DOMESTIC HOT WATER ***/
          ifVwater_tank: if (Vwater_tank > 0.0 ) then
                ! // convective heat flux to internal wall of hot water tank
                Qconv_water_to_inttankwall = &
                indoorConvectionHeatTransfer &
                (conv_coeff_intwall_tank, Asurf_tank, Tintwall_tank, Twater_tank)

                ! // heat flux by conduction through wall of hot water tank
                Qcond_tankwall =             &
                wallConduction               &
                (conductivity_wall_tank, Asurf_tank, Tintwall_tank, Textwall_tank, thickness_tankwall)

                ! // convective heat flux for external wall of hot water tank
                Qconv_exttankwall_to_indair = &
                outdoorConvectionHeatTransfer &
                (conv_coeff_extwall_tank, Asurf_tank, Textwall_tank, Tair_ind)
                ! // radiative heat flux for external wall of hot water tank
                ! //TODO: Should expand to consider windows and floor as well.

                Qlw_net_exttankwall_to_intwallroof = &
                outdoorRadiativeHeatTransfer         &
                (BVF_tank, Asurf_tank, emissivity_extwall_tank, Textwall_tank, Tintwallroof) ! // to building walls
                Qlw_net_exttankwall_to_indoormass =  &
                outdoorRadiativeHeatTransfer         &
                (MVF_tank, Asurf_tank, emissivity_extwall_tank, Textwall_tank, Tindoormass) ! // to internal mass

                ! // heat input into water of hot water tank
                qhwt_timestep =                      &
                heating                              &
                (setTwater_tank, Twater_tank, heating_efficiency_water, maxheatingpower_water)

                ! //Heat release from hot water heating due to efficiency losses
                Qloss_efficiency_heating_water =     &
                additionalSystemHeatingEnergy(qhwt_timestep, heating_efficiency_water)
!
          endif ifVwater_tank
!
!
!
!
          ifVwater_vessel: if (Vwater_vessel > 0.0 ) then
!
                ! // heat flux to internal wall of vessels holding DHW in use in building
                Qconv_water_to_intvesselwall =  &
                indoorConvectionHeatTransfer    &
                (conv_coeff_intwall_vessel, Awater_vessel, Tintwall_vessel, Twater_vessel)

                ! // heat flux by conduction through wall of vessels holding DHW in use in building
                Qcond_vesselwall =              &
                wallConduction                  &
                (conductivity_wall_vessel, Awater_vessel, Tintwall_vessel, Textwall_vessel, thickness_wall_vessel)

                ! // convective heat flux to external wall of vessels holding DHW in use in building
                Qconv_extvesselwall_to_indair = &
                outdoorConvectionHeatTransfer   &
                (conv_coeff_extwall_vessel, Awater_vessel, Textwall_vessel, Tair_ind)
                ! // radiative heat flux to external wall of vessels holding DHW in use in building
                ! // TODO: Should expand to consider windows and floor as well.
                Qlw_net_extvesselwall_to_wallroof = &
                outdoorRadiativeHeatTransfer        &
                (BVF_tank, Awater_vessel, emissivity_extwall_vessel, Textwall_vessel, Tintwallroof)
                Qlw_net_extvesselwall_to_indoormass = &
                outdoorRadiativeHeatTransfer          &
                (MVF_tank, Awater_vessel, emissivity_extwall_vessel, Textwall_vessel, Tindoormass)

                ! //Heat transfer due to use and replacement of water
                ! //qhwt_v = waterUseHeatTransfer(density_water, cp_water, flowrate_water_supply, Tincomingwater_tank, Twater_tank)
!
          elseif ( Vwater_vessel == minVwater_vessel .and.              &
           flowrate_water_supply < flowrate_water_drain ) then
!
                ! //Set Drain Flow rate to be same as usage flow rate to avoid hot water storage going below the set minimum threshold (minVwater_vessel)
                flowrate_water_drain = flowrate_water_supply
!
          endif ifVwater_vessel
!
!
!
!
!            //TODO: Need to consider the latent component of heat transfer for the DHW in use here also.

!            /************************************************************/
!            // Need to calculate the temperature change in the hot water tank due to energy flux
            Qtotal_net_water_tank = qhwt_timestep - Qconv_water_to_inttankwall
!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change in the hot water tank internal walls to energy flux
            Qtotal_net_intwall_tank = Qconv_water_to_inttankwall - Qcond_tankwall
!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change in the hot water tank external walls to energy flux
            Qtotal_net_extwall_tank =                                               & 
            Qcond_tankwall - Qconv_exttankwall_to_indair -                          &
            Qlw_net_exttankwall_to_intwallroof - Qlw_net_exttankwall_to_indoormass
!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change in the DHW vessels (i.e. in use hot water) due to energy flux
            Qtotal_net_water_vessel = - Qconv_water_to_intvesselwall
!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change in tthe DHW vessels (i.e. in use hot water) internal walls to energy flux
            Qtotal_net_intwall_vessel = Qconv_water_to_intvesselwall - Qcond_vesselwall
!            /************************************************************/

!            /************************************************************/
!            // Need to calculate the temperature change in the hot water tank external walls to energy flux
            Qtotal_net_extwall_vessel =                                              &
            Qcond_vesselwall - Qconv_extvesselwall_to_indair -                       &
            Qlw_net_extvesselwall_to_wallroof - Qlw_net_extvesselwall_to_indoormass
!            /************************************************************/

!            /*** END DOMESTIC HOT WATER ***/
!            /******************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change in the internal mass object due to energy flux, given its heat capacity
            Qtotal_net_indoormass =                                                  &
            Qsw_transmitted_window + Qconv_indair_to_indoormass +                    &
            Qlw_net_intwallroof_to_allotherindoorsurfaces +                          &
            Qlw_net_exttankwall_to_indoormass + Qlw_net_extvesselwall_to_indoormass
!
!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change in the internal volume of air due to energy flux to/from internal mass object AND to/from internal walls.
!            // and to/from hot water tanks and to/from the DHW held in use in the building.
            Qtotal_net_indair =                                                           &
            Q_appliance + Qmetabolic_sensible + Q_ventilation +                           &
            q_heating_timestep - q_cooling_timestep - Qconv_indair_to_indoormass -        &
            Qlw_net_intwallroof_to_allotherindoorsurfaces - Qconv_indair_to_intwallroof - &
            Qconv_indair_to_intwindow - Qconv_indair_to_intgroundfloor +                  &
            Qloss_efficiency_heating_air + Qconv_exttankwall_to_indair +                  &
            Qconv_extvesselwall_to_indair + Qloss_efficiency_heating_water
!
!            // TODO: Have not yet considered the latent heat component from occupancy gains!!
!
!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change of the internal wall surface due to internal heat exchanges and conduction through wall. Need an internal wall thermal mass also.
            Qtotal_net_intwallroof =                                                &
            Qconv_indair_to_intwallroof - Qcond_wallroof -                          &
            Qlw_net_intwallroof_to_allotherindoorsurfaces +                         &
            Qlw_net_exttankwall_to_intwallroof + Qlw_net_extvesselwall_to_wallroof
!
!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change of the external wall surface due to the conduction through the wall as well as wall surface exchanges with outside environment. Need an external wall thermal mass also.
            Qtotal_net_extwallroof =                                     &
            Qcond_wallroof + Qsw_absorbed_wallroof -                     &
            Qlw_net_extwallroof_to_outair - Qconv_extwallroof_to_outair

!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change of the internal window surface due to internal heat exchanges and conduction through wall. Need an internal window thermal mass also.
            Qtotal_net_intwindow =                       &
            Qconv_indair_to_intwindow - Qcond_window -   &
            Qlw_net_intwindow_to_allotherindoorsurfaces

!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change of the external window surface due to the conduction through the window as well as window surface exchanges with outside environment. Need an external window thermal mass also.
            Qtotal_net_extwindow =                                   &  
            Qcond_window + Qsw_absorbed_window -                     &
            Qlw_net_extwindow_to_outair - Qconv_extwindow_to_outair
!
!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change of the internal ground floor surface due to internal heat exchanges and conduction through floor. Need an internal floor thermal mass also.
            Qtotal_net_intgroundfloor =                           &
            Qconv_indair_to_intgroundfloor - Qcond_groundfloor -  &
            Qlw_net_intgroundfloor_to_allotherindoorsurfaces

!            /************************************************************/
!
!            /************************************************************/
!            // Need to calculate the temperature change of the external ground floor surface due to the conduction through the floor as well as surface exchanges with outside ground environment. Need an external wall thermal mass also.
            Qtotal_net_extgroundfloor = Qcond_groundfloor - Qcond_ground

!            // Accumulate the heat storage flux by wall/roof, window, floor, internal mass and air
            QS_total =                                               &
            Qtotal_net_extwallroof + Qtotal_net_intwallroof +        &
            Qtotal_net_extwindow + Qtotal_net_intwindow +            &
            Qtotal_net_extgroundfloor + Qtotal_net_intgroundfloor +  &
            Qtotal_net_indoormass + Qtotal_net_indair
!
            QS_fabric =                                                                     &
            Qtotal_net_extwallroof + Qtotal_net_intwallroof + Qtotal_net_extwindow +        &
            Qtotal_net_intwindow + Qtotal_net_extgroundfloor + Qtotal_net_intgroundfloor +  &
            Qtotal_net_indoormass
!
            QS_air = Qtotal_net_indair
!            /************************************************************/
!
!            /************************************************************/
!            // Detected Qf for SUEWS will be the outside wall and window heat flux, ventilation losses, and any additional HVAC system heat ejection related to system COP/efficiency.
!            // TODO: Consider Heat Storage Flux in this calculation. Where these fluxes are negative they represent the building absorbing heat from outside. This is part of Storage heat flux. There is also a question of how (or even if you can)  separate storage flux from QfB for solar passive gains.
!            //if (Qlw_net_extwallroof_to_outair < 0) Qlw_net_extwallroof_to_outair=0;  These lines are deleted, negatiev fluxes could happen, indicating building absorbs heat from outdoor environment
!            //if (Qlw_net_extwindow_to_outair < 0) Qlw_net_extwindow_to_outair=0;
!            //if (qwoc < 0) qwoc=0;
!            //if (Qconv_extwindow_to_outair < 0) Qconv_extwindow_to_outair=0;
!            //if (Q_ventilation > 0) Q_ventilation=0;
!
            Qf_ground_timestep = Qcond_ground * resolution
!
!            /************************************************************/
!
!            /************************************************************/
!            // STS - 31/0718: Adds time resolution heat flux to all building Energy Exchanges
!            // for building model timestep. Allows for calculating QfB and
!            // Qs outside of building energy model.
            Qsw_transmitted_window_tstepTotal =                                      &
            Qsw_transmitted_window_tstepTotal + Qsw_transmitted_window * resolution
!
            Qsw_absorbed_window_tstepTotal =                                         &
            Qsw_absorbed_window_tstepTotal + Qsw_absorbed_window * resolution
!
            Qsw_absorbed_wallroof_tstepTotal =                                       &
            Qsw_absorbed_wallroof_tstepTotal + Qsw_absorbed_wallroof * resolution
!
            Qconv_indair_to_indoormass_tstepTotal =                                         &
            Qconv_indair_to_indoormass_tstepTotal + Qconv_indair_to_indoormass * resolution
!
            Qlw_net_intwallroof_to_allotherindoorsurfaces_tstepTotal =  &
            Qlw_net_intwallroof_to_allotherindoorsurfaces_tstepTotal +  &
            Qlw_net_intwallroof_to_allotherindoorsurfaces * resolution
!
            Qlw_net_intwindow_to_allotherindoorsurfaces_tstepTotal =    &
            Qlw_net_intwindow_to_allotherindoorsurfaces_tstepTotal +    &
            Qlw_net_intwindow_to_allotherindoorsurfaces * resolution
!
            Qlw_net_intgroundfloor_to_allotherindoorsurfaces_tstepTotal =  &
            Qlw_net_intgroundfloor_to_allotherindoorsurfaces_tstepTotal +  &
            Qlw_net_intgroundfloor_to_allotherindoorsurfaces * resolution
!
            Q_appliance_tstepTotal = Q_appliance_tstepTotal + Q_appliance * resolution
!
            Q_ventilation_tstepTotal = Q_ventilation_tstepTotal + Q_ventilation * resolution
!
            Qconv_indair_to_intwallroof_tstepTotal =  &
            Qconv_indair_to_intwallroof_tstepTotal +  &
            Qconv_indair_to_intwallroof * resolution
!
            Qconv_indair_to_intwindow_tstepTotal =                                         &
            Qconv_indair_to_intwindow_tstepTotal + Qconv_indair_to_intwindow * resolution
!
            Qconv_indair_to_intgroundfloor_tstepTotal =                                              &
            Qconv_indair_to_intgroundfloor_tstepTotal + Qconv_indair_to_intgroundfloor * resolution

            Qloss_efficiency_heating_air_tstepTotal =                                            &
            Qloss_efficiency_heating_air_tstepTotal + Qloss_efficiency_heating_air * resolution

            Qcond_wallroof_tstepTotal = Qcond_wallroof_tstepTotal + Qcond_wallroof * resolution

            Qcond_window_tstepTotal = Qcond_window_tstepTotal + Qcond_window * resolution

            Qcond_groundfloor_tstepTotal = Qcond_groundfloor_tstepTotal + Qcond_groundfloor * resolution

            Qcond_ground_tstepTotal = Qcond_ground_tstepTotal + Qcond_ground * resolution

            Qlw_net_extwallroof_to_outair_tstepTotal =                                             &
            Qlw_net_extwallroof_to_outair_tstepTotal + Qlw_net_extwallroof_to_outair * resolution

            Qlw_net_extwindow_to_outair_tstepTotal =                                           &
            Qlw_net_extwindow_to_outair_tstepTotal + Qlw_net_extwindow_to_outair * resolution

            Qconv_extwallroof_to_outair_tstepTotal =                                           &
            Qconv_extwallroof_to_outair_tstepTotal + Qconv_extwallroof_to_outair * resolution

            Qconv_extwindow_to_outair_tstepTotal =                                         &
            Qconv_extwindow_to_outair_tstepTotal + Qconv_extwindow_to_outair * resolution

            q_cooling_timestepTotal =                                                      &
            q_cooling_timestepTotal +                                                      &
            (q_cooling_timestep + additionalSystemCoolingEnergy(q_cooling_timestep, coeff_performance_cooling)) * resolution

            qsensible_timestepTotal = qsensible_timestepTotal + Qmetabolic_sensible*resolution

            qlatent_timestepTotal = qlatent_timestepTotal + Qmetabolic_latent*resolution
!            /************************************************************/
!            //Adds timestep of heat storage flux
            QS_tstepTotal = QS_tstepTotal + QS_total * resolution
            QS_fabric_tstepTotal = QS_fabric_tstepTotal + QS_fabric*resolution
            QS_air_tstepTotal = QS_air_tstepTotal + QS_air*resolution
!            /************************************************************/
!            // Adds timestep heating/cooling to overall heating/cooling for building model
            Qtotal_heating = Qtotal_heating + (q_heating_timestep * resolution)
            Qtotal_cooling = Qtotal_cooling + (q_cooling_timestep * resolution)
!
!            /***** DOMESTIC HOT WATER HEATING *****/
            Qtotal_water_tank = Qtotal_water_tank + (qhwt_timestep * resolution)
!            /**************************************/
!
!            /************************************************************/
!
!            /*******************************************************************/
!            //Building temperature changes calculated for each lumped component
!            /*******************************************************************/
!            /************************************************************/
!            // Temperature Changes calculated for next timestep
!
!            /****************************/
!            /**** DOMESTIC HOT WATER ****/
!
!            // temperature (K) of DHW in use in Building due to heat transfer to building
            if ( Vwater_vessel > 0.0 ) then
                dTwater_vessel =                                                                      &
                (Qtotal_net_water_vessel / (density_water * cp_water) * Vwater_vessel) * resolution

                Twater_vessel = Twater_vessel + dTwater_vessel
            endif
!
!            //DHW "vessel" Internal Wall surface temperature (K)
            if (Vwall_vessel > 0.0) then
                dTintwall_vessel = &
                (Qtotal_net_intwall_vessel/((density_wall_vessel * cp_wall_vessel) * (Vwall_vessel / 2))) * resolution
                Tintwall_vessel = Tintwall_vessel + dTintwall_vessel

!            //DHW "vessel" External Wall surface temperature (K)
                dTextwall_vessel = &
                (Qtotal_net_extwall_vessel/((density_wall_vessel * cp_wall_vessel) * (Vwall_vessel / 2))) * resolution
                Textwall_vessel = Textwall_vessel + dTextwall_vessel
            endif
!
!            //Heat transfer to sewer/drain based on water in building going to drain.
            Qloss_drain =  &
            waterUseEnergyLossToDrains(density_water, cp_water, flowrate_water_drain, Twater_vessel, resolution)
!
!            //Need to recalculate the volume of DHW in use in building before calculating the temperature (K) of DHW in use.
            dVwater_vessel = (flowrate_water_supply - flowrate_water_drain) * resolution
            Vwater_vessel = Vwater_vessel + dVwater_vessel
!
!            //Avoid going below a given minimum volume
            if (Vwater_vessel < minVwater_vessel) then
                Vwater_vessel = minVwater_vessel
            endif
!
!            //Need to recalculate the temperature after mixing of water
            Twater_vessel =                                                                             &
            (((flowrate_water_supply * resolution) * Twater_tank) +                                     &
             ((Vwater_vessel - (flowrate_water_supply * resolution)) * Twater_vessel)) / Vwater_vessel
!
!            /**********/
!            //Recalculate DHW vessel surface area and wall volume based on new volume and maintaining the previous
            if (Vwater_vessel > 0.0)then ! //Checks that Volume isn't zero
                Awater_vessel = Vwater_vessel / VARatio_water_vessel
            else ! // if Volume is zero it makes the area also zero
                Awater_vessel = 0.0
            endif
            Vwall_vessel = Awater_vessel * thickness_wall_vessel
!            /**********/
!
!            //Hot Water Tank water temperature (K)
            if (Vwater_tank > 0.0) then
                dTwater_tank = (Qtotal_net_water_tank/((density_water * cp_water) * Vwater_tank)) * resolution !//(Q_hwt/((density_water*cp_water)*Vwater_tank))*resolution ! // resolution in seconds
                Twater_tank = Twater_tank + dTwater_tank
!//                write(*,*)"HWT Water Temperature: %f, dTwt: %f, Heating Q_hwt_timestep: %f \n", Twater_tank, dTwt, qhwt_timestep
            endif
!
!            // Hot Water Tank Internal Wall surface temperature (K)
            dTintwall_tank =                                                                                &
            (Qtotal_net_intwall_tank/((density_wall_tank * cp_wall_tank) * (Vwall_tank / 2))) * resolution
            Tintwall_tank = Tintwall_tank + dTintwall_tank
!
!            // Hot Water Tank External Wall surface temperature (K)
            dTextwall_tank =                                                                                &
            (Qtotal_net_extwall_tank/((density_wall_tank * cp_wall_tank) * (Vwall_tank / 2))) * resolution
            Textwall_tank = Textwall_tank + dTextwall_tank
!
!            //Need to recalculate the volume of water in hot water tank before calculating the temperature (K) after mixing with mains water
!            //NOTE: Currently not implemented as the water tank volume is considered constant at all times.
!            //dVwt = (flowrate_water_supply - flowrate_water_supply)*resolution;
!            //Vwater_tank = Vwater_tank + dVwt;
!
!            //Need to recalculate the temperature after mixing of mains water with hot water remaining in the tank
            Twater_tank =                                                                         &
            (((flowrate_water_supply * resolution) * Tincomingwater_tank) +                       &
             ((Vwater_tank - (flowrate_water_supply * resolution)) * Twater_tank)) / Vwater_tank
!
!            /** END DOMESTIC HOT WATER **/
!            /****************************/
!
!
!            // Indoor thermal mass temperature (K)
            dTindoormass = (Qtotal_net_indoormass / ((density_indoormass * cp_indoormass) * Vindoormass)) * resolution ! // resolution in seconds
            Tindoormass = Tindoormass + dTindoormass
!
!            // Indoor air temperature (K)
            dTair_ind = (Qtotal_net_indair / ((density_air_ind * cp_air_ind) * Vair_ind)) * resolution ! // resolution in seconds
            Tair_ind = Tair_ind + dTair_ind
!
!            // print "dTi: " + str(dTi)
!            // print "Tair_ind: " + str(self.Tair_ind)
!            // print "Qtotal_net_indair: " + str(Qtotal_net_indair)
!            // print "q_heating_timestep: " + str(q_heating_timestep)
!
!            // Internal wall surface temperature (K), use x1 to split heat capacity, impacting surface temperature change with time
            dTintwallroof =                                                             & 
            (Qtotal_net_intwallroof / ((density_wallroof * cp_wallroof) *               &
             (Vwallroof * (1 - weighting_factor_heatcapacity_wallroof)))) * resolution  ! // resolution in seconds
            Tintwallroof = Tintwallroof + dTintwallroof
!
!            // print "Tintwallroof: " + str(self.Tintwallroof)
!
!            // External wall surface temperature (K)
            dTextwallroof =                                                       &
            (Qtotal_net_extwallroof / ((density_wallroof * cp_wallroof) *         &
             (Vwallroof * weighting_factor_heatcapacity_wallroof))) * resolution   !  // resolution in seconds
            Textwallroof = Textwallroof + dTextwallroof
!
!            // Internal window surface temperature (K)
            dTintwindow = (Qtotal_net_intwindow / ((density_window * cp_window) * (Vwindow / 2))) * resolution ! // resolution in seconds
            Tintwindow = Tintwindow + dTintwindow
!
!            // External window surface temperture
            dTextwindow = (Qtotal_net_extwindow / ((density_window * cp_window) * (Vwindow / 2))) * resolution ! // resolution in seconds
            Textwindow = Textwindow + dTextwindow
!
!            // Internal ground floor surface temperature (K)
            dTintgroundfloor =                                                                             &
            (Qtotal_net_intgroundfloor / ((density_groundfloor * cp_groundfloor) * (Vgroundfloor / 2))) *  &
            resolution ! // resolution in seconds
            Tintgroundfloor = Tintgroundfloor + dTintgroundfloor
!
!            // External ground floor surface temperature (K)
            dTextgroundfloor =                                                                             &
            (Qtotal_net_extgroundfloor / ((density_groundfloor * cp_groundfloor) * (Vgroundfloor / 2))) *  &
            resolution !  // resolution in seconds
            Textgroundfloor = Textgroundfloor + dTextgroundfloor
!            /************************************************************/

        enddo looptime
    else !iftimestepresolution
!        printf("Timestep: %i not equally divisible by given resolution: %i.\n", timestep, resolution)
    endif
!
!
!
  end subroutine tstep
!
!
!
!
  subroutine reinitialiseTemperatures
  end subroutine reinitialiseTemperatures
!
!
!
!
  subroutine readnml(fnml,self)
!
  use modulestebbsprecision
  use modulestebbs, only : LBM
!
  implicit none
!
  type(LBM)                     :: self
  character(len=256),intent(in) :: fnml
!
  character(len=256) :: BuildingType, BuildingName
!
  real(rprc) ::                                  &
  Height,                                        &
  FootprintArea,                                 &
  RatioInternalVolume,                           &
  GroundDepth,                                   &
!  ratio_window_wall,                             &
  WWR,                                           &
  WallThickness,                                 &
  FloorThickness,                                &
  WindowThickness,                               &
  WallInternalConvectionCoefficient,             &
  InternalMassConvectionCoefficient,             &
  FloorInternalConvectionCoefficient,            &
  WindowInternalConvectionCoefficient,           &
  WallExternalConvectionCoefficient,             &
  WindowExternalConvectionCoefficient,           &
  WallEffectiveConductivity,                     &
  GroundFloorEffectiveConductivity,              &
  WindowEffectiveConductivity,                   &
  ExternalGroundConductivity,                    &
  WallDensity,                                   &
  Wallx1,                                        &
  WallExternalArea,                              &
  GroundFloorDensity,                            &
  WindowDensity,                                 &
  InternalMassDensity,                           &
  IndoorAirDensity,                              &
  WallCp,                                        &
  GroundFloorCp,                                 &
  WindowCp,                                      &
  InternalMassCp,                                &
  IndoorAirCp,                                   &
  WallExternalEmissivity,                        &
  WallInternalEmissivity,                        &
  InternalMassEmissivity,                        &
  WindowExternalEmissivity,                      &
  WindowInternalEmissivity,                      &
  WindowTransmissivity,                          &
  WindowAbsorbtivity,                            &
  WindowReflectivity,                            &
  WallTransmissivity,                            &
  WallAbsorbtivity,                              &
  WallReflectivity,                              &
  WallBuildingViewFactor,                        &
  WallGroundViewFactor,                          &
  WallSkyViewFactor,                             &
  occupants,                                     &
  MetabolicRate,                                 &
  ApplianceRating,                               &
  LatentSensibleRatio,                           &
  appliance_power_rating,                        &
  TotalNumberofAppliances,                       &
  ApplianceUsageFactor,                          &
  MaxHeatingPower,                               &
  HeatingSystemEfficiency,                       &
  MaxCoolingPower,                               &
  CoolingSystemCOP,                              &
  AirVolume,                                     &
  VentilationRate,                               &
  WallArea,                                      &
  WallVolume,                                    &
  FloorArea,                                     &
  FloorVolume,                                   &
  WindowArea,                                    &
  WindowVolume,                                  &
  InternalMassVolume,                            &
  InternalMassArea,                              &
  IndoorAirStartTemperature,                     &
  IndoorMassStartTemperature,                    &
  WallIndoorSurfaceTemperature,                  &
  WallOutdoorSurfaceTemperature,                 &
  WindowIndoorSurfaceTemperature,                &
  WindowOutdoorSurfaceTemperature,               &
  GroundFloorIndoorSurfaceTemperature,           &
  GroundFloorOutdoorSurfaceTemperature,          &
  HeatingSetpointTemperature,                    &
  CoolingSetpointTemperature,                    &
  WaterTankTemperature,                          &
  InternalWallWaterTankTemperature,              &
  ExternalWallWaterTankTemperature,              &
  WaterTankWallThickness,                        &
  MainsWaterTemperature,                         &
  WaterTankWaterVolume,                          &
  WaterTankSurfaceArea,                          &
  HotWaterTankWallVolume,                        &
  HotWaterHeatingSetpointTemperature,            &
  DomesticHotWaterTemperatureInUseInBuilding,    &
  InternalWallDHWVesselTemperature,              &
  ExternalWallDHWVesselTemperature,              &
  DHWVesselWallThickness,                        &
  DHWWaterVolume,                                &
  DHWSurfaceArea,                                &
  DHWVesselWallVolume,                           &
  DHWVesselEmissivity,                           &
  HotWaterFlowRate,                              &
  DHWDrainFlowRate,                              &
  DHWSpecificHeatCapacity,                       &
  HotWaterTankSpecificHeatCapacity,              &
  DHWVesselSpecificHeatCapacity,                 &
  DHWDensity,                                    &
  HotWaterTankWallDensity,                       &
  DHWVesselDensity,                              &
  HotWaterTankBuildingWallViewFactor,            &
  HotWaterTankInternalMassViewFactor,            &
  HotWaterTankWallConductivity,                  &
  HotWaterTankInternalWallConvectionCoefficient, &
  HotWaterTankExternalWallConvectionCoefficient, &
  HotWaterTankWallEmissivity,                    &
  DHWVesselWallConductivity,                     &
  DHWVesselInternalWallConvectionCoefficient,    &
  DHWVesselExternalWallConvectionCoefficient,    &
  DHWVesselWallEmissivity,                       &
  MaximumHotWaterHeatingPower,                   &
  HotWaterHeatingEfficiency,                     &
  MinimumVolumeOfDHWinUse
!
  namelist/specification/                     &
  BuildingType,                               &
  BuildingName,                               &
  Height,                                     &
  FootprintArea,                              &
  WallExternalArea,                           &
  WWR,                                        &
  RatioInternalVolume,                        &
  WallThickness,                              &
  WallEffectiveConductivity,                  &
  WallDensity,                                &
  WallCp,                                     &
  Wallx1,                                     &
  FloorThickness,                             &
  GroundFloorEffectiveConductivity,           &
  GroundFloorDensity,                         &
  GroundFloorCp,                              &
  GroundDepth,                                &
  ExternalGroundConductivity,                 &
  WindowThickness,                            &
  WindowEffectiveConductivity,                &
  WindowDensity,                              &
  WindowCp,                                   &
  InternalMassDensity,                        &
  InternalMassCp,                             &
  IndoorAirDensity,                           &
  IndoorAirCp,                                &
  WallInternalConvectionCoefficient,          &
  InternalMassConvectionCoefficient,          &
  FloorInternalConvectionCoefficient,         &
  WindowInternalConvectionCoefficient,        &
  WallExternalConvectionCoefficient,          &
  WindowExternalConvectionCoefficient,        &
  WallExternalEmissivity,                     &
  WallInternalEmissivity,                     &
  InternalMassEmissivity,                     &
  WindowExternalEmissivity,                   &
  WindowInternalEmissivity,                   &
  WindowTransmissivity,                       &
  WindowAbsorbtivity,                         &
  WindowReflectivity,                         &
  WallTransmissivity,                         &
  WallAbsorbtivity,                           &
  WallReflectivity,                           &
  WallBuildingViewFactor,                     &
  WallGroundViewFactor,                       &
  WallSkyViewFactor,                          &
  Occupants,                                  &
  MetabolicRate,                              &
  LatentSensibleRatio,                        &
  ApplianceRating,                            &
  TotalNumberofAppliances,                    &
  ApplianceUsageFactor,                       &
  MaxHeatingPower,                            &
  HeatingSystemEfficiency,                    &
  MaxCoolingPower,                            &
  CoolingSystemCOP,                           &
  VentilationRate,                            &
  IndoorAirStartTemperature,                  &
  IndoorMassStartTemperature,                 &
  WallIndoorSurfaceTemperature,               &
  WallOutdoorSurfaceTemperature,              &
  WindowIndoorSurfaceTemperature,             &
  WindowOutdoorSurfaceTemperature,            &
  GroundFloorIndoorSurfaceTemperature,        &
  GroundFloorOutdoorSurfaceTemperature,       &
  HeatingSetpointTemperature,                 &
  CoolingSetpointTemperature,                 &
  WaterTankTemperature,                       &
  InternalWallWaterTankTemperature,           &
  ExternalWallWaterTankTemperature,           &
  WaterTankWallThickness,                     &
  MainsWaterTemperature,                      &
  WaterTankWaterVolume,                       &
  WaterTankSurfaceArea,                       &
  HotWaterHeatingSetpointTemperature,         &
  HotWaterTankWallEmissivity,                 &
  DomesticHotWaterTemperatureInUseInBuilding, &
  InternalWallDHWVesselTemperature,           &
  ExternalWallDHWVesselTemperature,           &
  DHWVesselWallThickness,                     &
  DHWWaterVolume,                             &
  DHWSurfaceArea,                             &
  DHWVesselEmissivity,                        &
  HotWaterFlowRate,                           &
  DHWDrainFlowRate,                           &
  DHWSpecificHeatCapacity,                    &
  HotWaterTankSpecificHeatCapacity,           &
  DHWVesselSpecificHeatCapacity,              &
  DHWDensity,                                 &
  HotWaterTankWallDensity,                    &
  DHWVesselDensity,                           &
  HotWaterTankBuildingWallViewFactor,         &
  HotWaterTankInternalMassViewFactor,         &
  HotWaterTankWallConductivity,               &
  HotWaterTankInternalWallConvectionCoefficient, &
  HotWaterTankExternalWallConvectionCoefficient, &
  DHWVesselWallConductivity,                     &
  DHWVesselInternalWallConvectionCoefficient,    &
  DHWVesselExternalWallConvectionCoefficient,    &
  DHWVesselWallEmissivity,                       &
  MaximumHotWaterHeatingPower,                   &
  HotWaterHeatingEfficiency,                     &
  MinimumVolumeOfDHWinUse
!
!
!
! Maybe namelist nml can be more general
!
  open(8,file=trim(fnml) )
  read(8,nml=specification)
!
!
!
! Asign to the object
!
  self%BuildingType        = BuildingType
  self%BuildingName        = BuildingName
  self%ratio_window_wall   = WWR
  self%Afootprint          = FootprintArea  ! # NEW
  self%height_building     = Height
  self%wallExternalArea    = WallExternalArea ! # NEW
  self%ratioInternalVolume = RatioInternalVolume ! # NEW ratio internal/external volume
  self%thickness_wallroof  = WallThickness ! # wall and roof thickness as not separate (in metres)
  self%thickness_groundfloor  = FloorThickness ! # ground floor thickness (in metres)
  self%depth_ground           = GroundDepth ! # depth of ground considered for QfB to ground (in metres)
  self%thickness_window       = WindowThickness ! # all window's thickness (in metres)
  self%conv_coeff_intwallroof = WallInternalConvectionCoefficient
  self%conv_coeff_indoormass  = InternalMassConvectionCoefficient
  self%conv_coeff_intgroundfloor = FloorInternalConvectionCoefficient
  self%conv_coeff_intwindow   = WindowInternalConvectionCoefficient
  self%conv_coeff_extwallroof = WallExternalConvectionCoefficient
  self%conv_coeff_extwindow   = WindowExternalConvectionCoefficient
  self%conductivity_wallroof  = WallEffectiveConductivity
  self%conductivity_groundfloor = GroundFloorEffectiveConductivity
  self%conductivity_window    = WindowEffectiveConductivity
  self%conductivity_ground    = ExternalGroundConductivity
  self%density_wallroof       = WallDensity
  self%weighting_factor_heatcapacity_wallroof = Wallx1 ! #to split heat capacity of external and internal node of wall/roof, impacting surface temperature change with time
  self%density_groundfloor    = GroundFloorDensity
  self%density_window         = WindowDensity
  self%density_indoormass     = InternalMassDensity
  self%density_air_ind        = IndoorAirDensity
  self%cp_wallroof            = WallCp
  self%cp_groundfloor         = GroundFloorCp
  self%cp_window              = WindowCp
  self%cp_indoormass          = InternalMassCp
  self%cp_air_ind             = IndoorAirCp
  self%emissivity_extwallroof = WallExternalEmissivity
  self%emissivity_intwallroof = WallInternalEmissivity
  self%emissivity_indoormass  = InternalMassEmissivity
  self%emissivity_extwindow   = WindowExternalEmissivity
  self%emissivity_intwindow   = WindowInternalEmissivity
  self%windowTransmissivity   = WindowTransmissivity
  self%windowAbsorbtivity     = WindowAbsorbtivity
  self%windowReflectivity     = WindowReflectivity
  self%wallTransmisivity      = WallTransmissivity
  self%wallAbsorbtivity       = WallAbsorbtivity
  self%wallReflectivity       = WallReflectivity
  self%BVF_extwall            = WallBuildingViewFactor
  self%GVF_extwall            = WallGroundViewFactor
  self%SVF_extwall            = WallSkyViewFactor
  self%occupants              = occupants ! #0
!  # self.maxOccupants = nameListFile["Occupants"]
!  # self.appliance_energy = 0
!  # self.metabolic_energy = 0
!  # self.water_users = 0
!  # self.BuildingClass = nameListFile["BuildingClass"]
!  # self.BuildingCount = nameListFile["BuildingCount"]
   self%metabolic_rate        = MetabolicRate
   self%ratio_metabolic_latent_sensible = LatentSensibleRatio
   self%appliance_power_rating = ApplianceRating
   self%appliance_totalnumber  = int(TotalNumberofAppliances)
   self%appliance_usage_factor = ApplianceUsageFactor
   self%maxheatingpower_air    = MaxHeatingPower
   self%heating_efficiency_air = HeatingSystemEfficiency
   self%maxcoolingpower_air    = MaxCoolingPower
   self%coeff_performance_cooling = CoolingSystemCOP
   self%Vair_ind = (self%Afootprint * self%height_building) * &
                   (1 - self%ratioInternalVolume) ! # Multiplied by factor that accounts for internal mass
!   self%ventilation_rate       = VentilationRate ! # Fixed at begining to have no natural ventilation. Given in units of volume of air per second
   self%ventilation_rate       = self%Vair_ind * VentilationRate / 3600.0 ! Now Nakao includes setVentilationRate(self, VR: float) in STEBBS.py
   self%Awallroof              = (self%wallExternalArea * (1 - self%ratio_window_wall)) + self%Afootprint ! # last component accounts for the roof as not considered seperately in the model
   self%Vwallroof              = self%Awallroof * self%thickness_wallroof
   self%Vgroundfloor           = self%Afootprint * self%thickness_groundfloor
   self%Awindow                = self%wallExternalArea * self%ratio_window_wall
   self%Vwindow                = self%Awindow * self%thickness_window
   self%Vindoormass            = self%Vair_ind * self%ratioInternalVolume ! # Multiplied by factor that accounts for internal mass as proportion of total air volume
   self%Aindoormass = 6 * (self%Vindoormass ** (2. / 3.)) ! # Assumed internal mass as a cube
   self%h_i   = (/self%conv_coeff_intwallroof, self%conv_coeff_indoormass, &
                  self%conv_coeff_intgroundfloor, self%conv_coeff_intwindow/)
   self%h_o   = (/self%conv_coeff_extwallroof, self%conv_coeff_extwindow/)
   self%k_eff = (/self%conductivity_wallroof, self%conductivity_groundfloor, &
                  self%conductivity_window, self%conductivity_ground/)
   self%rho   = (/self%density_wallroof, self%density_groundfloor, &
                  self%density_window, self%density_indoormass,    &
                  self%density_air_ind/)
   self%Cp    = (/self%cp_wallroof, self%cp_groundfloor, self%cp_window, &
                  self%cp_indoormass, self%cp_air_ind/)
   self%emis  = (/self%emissivity_extwallroof, self%emissivity_intwallroof, &
                  self%emissivity_indoormass, self%emissivity_extwindow,    &
                  self%emissivity_intwindow/)
   self%wiTAR = (/self%windowTransmissivity, self%windowAbsorbtivity, self%windowReflectivity/)
   self%waTAR = (/self%wallTransmisivity, self%wallAbsorbtivity, self%wallReflectivity/)

   self%viewFactors  = (/self%BVF_extwall, self%GVF_extwall, self%SVF_extwall/) ! # Building, ground, and sky view factors
   self%occupantData = (/self%occupants, self%metabolic_rate, self%ratio_metabolic_latent_sensible/)

!# List of occupant related factors [Number of occupants, average metabolic rate, latent to sensible heat ratio]
! #            self.applianceData = [self.appliance_power_rating,self.appliance_totalnumber,self.appliance_usage_factor] # List of appliance related factors [Average appliance power rating, Number of appliances, factor usage of appliances (0 to 1)]

! #            self.heatingSystem = [self.maxheatingpower_air,self.heating_efficiency_air]
! #            self.coolingSystem = [self.maxcoolingpower_air,self.coeff_performance_cooling]

  self%Tair_ind = IndoorAirStartTemperature + 273.15 ! # Indoor air temperature (K)
  self%Tindoormass = IndoorMassStartTemperature + 273.15 ! # Indoor mass temperature (K)
  self%Tintwallroof = WallIndoorSurfaceTemperature + 273.15 ! # Wall indoor surface temperature (K)
  self%Textwallroof = WallOutdoorSurfaceTemperature + 273.15 ! # Wall outdoor surface temperature (K)
  self%Tintwindow = WindowIndoorSurfaceTemperature + 273.15 ! # Window indoor surface temperature (K)
! # Window outdoor surface temperature (K)
  self%Textwindow = WindowOutdoorSurfaceTemperature + 273.15
! # Ground floor indoor surface temperature (K)
  self%Tintgroundfloor = GroundFloorIndoorSurfaceTemperature + 273.15
! # Ground floor outdoor surface temperature (K)
  self%Textgroundfloor = GroundFloorOutdoorSurfaceTemperature + 273.15
! # Heating and Cooling setpoint temperature (K)s, respectively

  self%Ts     = (/HeatingSetpointTemperature + 273.15,  &
                  CoolingSetpointTemperature + 273.15/)
!
  self%initTs = (/HeatingSetpointTemperature + 273.15,  &
                  CoolingSetpointTemperature + 273.15/)
!
  self%HTsAverage  = (/18 + 273.15, 18 + 273.15, 18 + 273.15/)
  self%HWTsAverage = (/10 + 273.15, 10 + 273.15, 10 + 273.15/)

!            ####################################################
!            ################ DOMESTIC HOT WATER ################

  self%Twater_tank = WaterTankTemperature + 273.15 ! # Water temperature (K) in Hot Water Tank
  self%Tintwall_tank = InternalWallWaterTankTemperature + 273.15 ! # Hot water tank internal wall temperature (K)
  self%Textwall_tank = ExternalWallWaterTankTemperature + 273.15 ! # Hot water tank external wall temperature (K)
  self%thickness_tankwall = WaterTankWallThickness ! # Hot water tank wall thickness
  self%Tincomingwater_tank = MainsWaterTemperature + 273.15 ! # Water temperature (K) of Water coming into the Water Tank
  self%Vwater_tank = WaterTankWaterVolume ! # Volume of Water in Hot Water Tank (m^3)
  self%Asurf_tank = WaterTankSurfaceArea ! # Surface Area of Hot Water Tank(m^2)
  self%Vwall_tank = self%Asurf_tank * self%thickness_tankwall ! # Wall volume of Hot Water Tank(m^2)
  self%setTwater_tank = HotWaterHeatingSetpointTemperature + 273.15 ! # Water Tank setpoint temperature (K)
  self%init_wtTs = HotWaterHeatingSetpointTemperature + 273.15 ! # Initial Water Tank setpoint temperature (K)
  self%Twater_vessel = DomesticHotWaterTemperatureInUseInBuilding + 273.15 ! # Water temperature (K) of water held in use in Building
  self%Tintwall_vessel = InternalWallDHWVesselTemperature + 273.15 ! # Hot water vessel internal wall temperature (K)
  self%Textwall_vessel = ExternalWallDHWVesselTemperature + 273.15 ! # Hot water tank external wall temperature (K)
  self%thickness_wall_vessel = DHWVesselWallThickness ! # DHW vessels wall thickness

  self%Vwater_vessel = DHWWaterVolume ! # Volume of water held in use in building
  self%Awater_vessel = DHWSurfaceArea ! # Surface Area of Hot Water in Vessels in Building
  self%Vwall_vessel = self%Awater_vessel * self%thickness_wall_vessel ! # Wall volume of Hot water Vessels in Building
  self%flowrate_water_supply = HotWaterFlowRate ! # Hot Water Flow Rate in m^3 / s
  self%flowrate_water_drain  = DHWDrainFlowRate ! # Draining of Domestic Hot Water held in building

  self%single_flowrate_water_supply = HotWaterFlowRate ! # Hot Water Flow Rate in m^3 s^-1 for a single HW unit
  self%single_flowrate_water_drain  = DHWDrainFlowRate ! # Draining of Domestic Hot Water held in building

  self%cp_water                     = DHWSpecificHeatCapacity ! # Specific Heat Capacity of Domestic Hot Water
  self%cp_wall_tank   = HotWaterTankSpecificHeatCapacity ! # Specific Heat Capacity of Hot Water Tank wall
  self%cp_wall_vessel = DHWVesselSpecificHeatCapacity ! # Specific Heat Capacity of Vessels containing DHW in use in Building

  self%density_water     = DHWDensity ! # Density of water
  self%density_wall_tank = HotWaterTankWallDensity ! # Density of hot water tank wall

  self%density_wall_vessel = DHWVesselDensity ! # Density of vessels containing DHW in use in buildings

  self%BVF_tank = HotWaterTankBuildingWallViewFactor ! # water tank - building wall view factor
  self%MVF_tank = HotWaterTankInternalMassViewFactor ! # water tank - building internal mass view factor

  self%conductivity_wall_tank  = HotWaterTankWallConductivity ! # Effective Wall conductivity of the Hot Water Tank
  self%conv_coeff_intwall_tank = HotWaterTankInternalWallConvectionCoefficient ! # Effective Internal Wall convection coefficient of the Hot Water Tank
  self%conv_coeff_extwall_tank = HotWaterTankExternalWallConvectionCoefficient ! # Effective External Wall convection coefficient of the Hot Water Tank

  self%emissivity_extwall_tank = HotWaterTankWallEmissivity ! # Effective External Wall emissivity of the Hot Water Tank
  self%conductivity_wall_vessel = DHWVesselWallConductivity ! # Effective Conductivity of vessels containing DHW in use in Building.
  self%conv_coeff_intwall_vessel = DHWVesselInternalWallConvectionCoefficient ! # Effective Internal Wall convection coefficient of the Vessels holding DHW in use in Building
  self%conv_coeff_extwall_vessel = DHWVesselExternalWallConvectionCoefficient ! # Effective Enternal Wall convection coefficient of the Vessels holding DHW in use in Building
  self%emissivity_extwall_vessel = DHWVesselWallEmissivity ! # Effective External Wall emissivity of hot water being used within building

! ## NOTE THAT LATENT HEAT FLUX RELATING TO HOT WATER USE CURRENTLY NOT IMPLEMENTED ##
!
  self%maxheatingpower_water    = MaximumHotWaterHeatingPower
  self%heating_efficiency_water = HotWaterHeatingEfficiency
  self%minVwater_vessel         = MinimumVolumeOfDHWinUse
!
  self%minHeatingPower_DHW = MaximumHotWaterHeatingPower
  self%HeatingPower_DHW    = MaximumHotWaterHeatingPower
!
  self%HWPowerAverage = (/30000, 30000, 30000/)
!
!
!
  return
!
  end subroutine readnml
!
!
!
!
  subroutine create_building(case,self,icase)
!
  use modulestebbs, only : LBM
!
  implicit none
!
  integer,intent(in) :: icase 
  character(len=256) :: case
  type(LBM) :: self
!
!
!
  self%idLBM = icase
  self%fnmlLBM = './BuildClasses/' // trim(case) // '.nml'
  self%case = trim(case)
  self%Qtotal_heating = 0.0  ! # currently only sensible but this needs to be  split into sensible and latent heat components
  self%Qtotal_cooling = 0.0  ! # currently only sensible but this needs to be  split into sensible and latent heat components

  self%Qmetabolic_sensible = 0.0  ! # Sensible heat flux from people in building
  self%Qmetabolic_latent = 0.0  ! # Latent heat flux from people in building

!        ############ DOMESTIC HOT WATER TOTAL HEATING ############
  self%Qtotal_water_tank = 0.0
  self%qhwtDrain = 0.0
!
  self%EnergyExchanges(:) = 0.0
!        ########## END DOMESTIC HOT WATER TOTAL HEATING ##########
!
  ifnonml: if ( trim(case) == 'none' ) then
            self%BuildingType = 'None'
            self%BuildingName = 'Default'
            self%ratio_window_wall = 0.4  ! # window to wall ratio
            self%Afootprint = 225.0
            self%height_building = 15.0
            self%wallExternalArea = 450.0
            self%ratioInternalVolume = 0.1
            self%thickness_wallroof = 0.25  ! # wall and roof thickness as not separate (in metres)
            self%thickness_groundfloor = 0.5  ! # ground floor thickness (in metres)
            self%depth_ground = 2.0  ! # depth of ground considered for calculating QfB to ground (in metres)
            self%thickness_window = 0.05  ! # all window's thickness (in metres)
            self%conv_coeff_intwallroof = 1
            self%conv_coeff_indoormass = 1
            self%conv_coeff_intgroundfloor = 1
            self%conv_coeff_intwindow = 1
            self%conv_coeff_extwallroof = 2 ! # Could be changed to react to outdoor wind speed from SUEWS 
            self%conv_coeff_extwindow = 2
            self%conductivity_wallroof = 0.2
            self%conductivity_groundfloor = 0.2
            self%conductivity_window = 0.5
            self%conductivity_ground = 1.0
            self%density_wallroof = 50
            self%weighting_factor_heatcapacity_wallroof = 0.5
            self%density_groundfloor = 1000
            self%density_window = 5
            self%density_indoormass = 250
            self%density_air_ind = 1.225
            self%cp_wallroof = 500
            self%cp_groundfloor = 1500
            self%cp_window = 840
            self%cp_indoormass = 1000
            self%cp_air_ind = 1005
            self%emissivity_extwallroof = 0.85
            self%emissivity_intwallroof = 0.85
            self%emissivity_indoormass = 0.85
            self%emissivity_extwindow = 0.85
            self%emissivity_intwindow = 0.85
            self%windowTransmissivity = 0.9
            self%windowAbsorbtivity = 0.05
            self%windowReflectivity = 0.05
            self%wallTransmisivity = 0.0
            self%wallAbsorbtivity = 0.5
            self%wallReflectivity = 0.5
            self%BVF_extwall = 0.4
            self%GVF_extwall = 0.2
            self%SVF_extwall = 0.4
            self%occupants = 15
            self%metabolic_rate = 250
            self%ratio_metabolic_latent_sensible = 0.8
            self%appliance_power_rating = 100
            self%appliance_totalnumber = 45
            self%appliance_usage_factor = 0.8
            self%maxheatingpower_air = 45000
            self%heating_efficiency_air = 0.9
            self%maxcoolingpower_air = 3000
            self%coeff_performance_cooling = 2.0
            self%Vair_ind =                             &
            (self%Afootprint * self%height_building) *  &
            (1 - self%ratioInternalVolume) ! # Multiplied by factor that accounts for internal mass
            self%ventilation_rate = 0  ! # Fixed at begining to have no natural ventilation. Given in units of volume of air per second
            self%Awallroof =                                          &
            (self%wallExternalArea * (1 - self%ratio_window_wall)) +  &
            self%Afootprint ! # last component accounts for the roof as not considered seperately in the model
            self%Vwallroof = self%Awallroof * self%thickness_wallroof
            self%Vgroundfloor = self%Afootprint * self%thickness_groundfloor
            self%Awindow = self%wallExternalArea * self%ratio_window_wall
            self%Vwindow = self%Awindow * self%thickness_window
            self%Vindoormass = self%Vair_ind * self%ratioInternalVolume ! # Multiplied by factor that accounts for internal mass as proportion of total air volume
            self%Aindoormass = 6 * (self%Vindoormass ** (2. / 3.)) ! # Assumed internal mass as a cube
            self%h_i = (/self%conv_coeff_intwallroof, self%conv_coeff_indoormass, &
                         self%conv_coeff_intgroundfloor, self%conv_coeff_intwindow/)

            self%h_o = (/self%conv_coeff_extwallroof, self%conv_coeff_extwindow/)
            self%k_eff = (/self%conductivity_wallroof, self%conductivity_groundfloor, &
                          self%conductivity_window, self%conductivity_ground/)

            self%rho = (/self%density_wallroof, self%density_groundfloor,  &
                         self%density_window, self%density_indoormass,     &
                         self%density_air_ind/)
            self%Cp = (/self%cp_wallroof, self%cp_groundfloor, self%cp_window, &
                        self%cp_indoormass, self%cp_air_ind/)
            self%emis = (/self%emissivity_extwallroof, self%emissivity_intwallroof, &
                          self%emissivity_indoormass, self%emissivity_extwindow,    &
                          self%emissivity_intwindow/)
            self%wiTAR = (/self%windowTransmissivity, self%windowAbsorbtivity, self%windowReflectivity/)
            self%waTAR = (/self%wallTransmisivity, self%wallAbsorbtivity, self%wallReflectivity/)

            self%viewFactors = (/self%BVF_extwall, self%GVF_extwall, self%SVF_extwall/) !  # Building, ground, and sky view factors
            self%occupantData = (/self%occupants, self%metabolic_rate, &
                                 self%ratio_metabolic_latent_sensible/)
!            #            self.applianceData = [self.appliance_power_rating,self.appliance_totalnumber,self.appliance_usage_factor] # List of appliance related factors [Average appliance power rating, Number of appliances, factor usage of appliances (0 to 1)]

!            #            self.heatingSystem = [self.maxheatingpower_air,self.heating_efficiency_air]
!            #            self.coolingSystem = [self.maxcoolingpower_air,self.coeff_performance_cooling]

            self%Tair_ind = 20 + 273.15 ! # Indoor air temperature (K)
            self%Tindoormass = 20 + 273.15 ! # Indoor mass temperature (K)
            self%Tintwallroof = 20 + 273.15 ! # Wall indoor surface temperature (K)
            self%Textwallroof = 20 + 273.15 ! # Wall outdoor surface temperature (K)
            self%Tintwindow = 20 + 273.15 ! # Window indoor surface temperature (K)
            self%Textwindow = 20 + 273.15 ! # Window outdoor surface temperature (K)
            self%Tintgroundfloor = 20 + 273.15 ! # Ground floor indoor surface temperature (K)
            self%Textgroundfloor = 20 + 273.15 ! # Ground floor outdoor surface temperature (K)
            self%Ts = (/18 + 273.15, 24 + 273.15/) ! # Heating and Cooling setpoint temperatures (K), respectively
            self%initTs = (/18 + 273.15, 24 + 273.15/)
            self%HTsAverage = (/18 + 273.15, 18 + 273.15, 18 + 273.15/) ! #
            self%HWTsAverage = (/10 + 273.15, 10 + 273.15, 10 + 273.15/)
            self%Twater_tank = 70 + 273.15 ! # Water temperature (K) in Hot Water Tank
            self%Tintwall_tank = 70 + 273.15 ! # Hot water tank internal wall temperature (K)
            self%Textwall_tank = 30 + 273.15 ! # Hot water tank external wall temperature (K)
            self%thickness_tankwall = 0.1 ! # Hot water tank wall thickness (m)
            self%Tincomingwater_tank = 10 + 273.15 ! # Water temperature (K) of Water coming into the Water Tank
            self%Vwater_tank = 2 ! # Volume of Water in Hot Water Tank (m^3)  h = 1.5, (2/(1.5*3.14))^0.5 = r =
            self%Asurf_tank = 8.8 ! # Surface Area of Hot Water Tank(m^2) - cylinder h= 1.5
            self%Vwall_tank = self%Asurf_tank * self%thickness_tankwall ! # Wall volume of Hot Water Tank(m^2)
            self%setTwater_tank = 60 + 273.15 ! # Water Tank setpoint temperature (K)
            self%init_wtTs = 60 + 273.15 ! # Initial Water Tank setpoint temperature (K)

            self%Twater_vessel = 35 + 273.15 ! # Water temperature (K) of water held in use in Building
            self%Tintwall_vessel = 35 + 273.15 ! # Hot water vessel internal wall temperature (K)
            self%Textwall_vessel = 25 + 273.15 ! # Hot water vessel external wall temperature (K)
            self%thickness_wall_vessel = 0.005 ! # DHW vessels wall thickness (m)

            self%Vwater_vessel = 0 ! # Volume of water held in use in building (m^3)
            self%Awater_vessel = 10 ! # Surface Area of Hot Water in Vessels in Building (m^2)
            self%Vwall_vessel = self%Awater_vessel * self%thickness_wall_vessel ! # Wall volume of Hot water Vessels in Building
            self%flowrate_water_supply = 0 ! # Hot Water Flow Rate in m^3 / s
            self%flowrate_water_drain = 0 ! # Draining of Domestic Hot Water held in building

            self%single_flowrate_water_supply = 0 ! # Hot Water Flow Rate in m^3 s^-1 for a single HW unit
            self%single_flowrate_water_drain = 0 ! # Draining of Domestic Hot Water held in building

            self%cp_water = 4180.1 ! # Specific Heat Capacity of Domestic Hot Water (J/kg K)
            self%cp_wall_tank = 1000 ! # Specific Heat Capacity of Hot Water Tank wall
            self%cp_wall_vessel = 1900 ! # Specific Heat Capacity of Vessels containing DHW in use in Building (value here is based on MDPE)

            self%density_water = 1000  ! # Density of water
            self%density_wall_tank = 50 ! # Density of hot water tank wall
            self%density_wall_vessel = 930 ! # Density of vessels containing DHW in use in buildings

            self%BVF_tank = 0.2 ! # water tank - building wall view factor
            self%MVF_tank = 0.8 ! # water tank - building internal mass view factor

            self%conductivity_wall_tank = 0.1 ! # Effective Wall conductivity of the Hot Water Tank (based on polyurethan foam given in https://www.lsta.lt/files/events/28_jarfelt.pdf and from https://www.sciencedirect.com/science/article/pii/S0360544214011189?via%3Dihub)
            self%conv_coeff_intwall_tank = 243 ! # Effective Internal Wall convection coefficient of the Hot Water Tank (W/m2 . K) given in http://orbit.dtu.dk/fedora/objects/orbit:77843/datastreams/file_2640258/content

            self%conv_coeff_extwall_vessel = 4 ! # Effective Enternal Wall convection coefficient of the Vessels holding DHW in use in Building
            self%emissivity_extwall_vessel = 0.88 ! # Effective External Wall emissivity of hot water being used within building

!            ## NOTE THAT LATENT HEAT FLUX RELATING TO HOT WATER USE CURRENTLY NOT IMPLEMENTED ##
!
            self%maxheatingpower_water = 3000 ! # Watts
            self%heating_efficiency_water = 0.95
            self%minVwater_vessel = 0.1 ! # m3

            self%minHeatingPower_DHW = 3000
            self%HeatingPower_DHW = 3000

            self%HWPowerAverage = (/30000, 30000, 30000/)


  else
          call readnml(self%fnmlLBM,self)
!
  endif ifnonml
!
!
!
  return
!
  end subroutine create_building
!
!
!
!
#ifndef STEBBSONLINE
  program main
#else
  subroutine main_for_offline
#endif
!
  use modulestebbs
!
  implicit none
!
  integer      :: i, j
!
!
!
  if ( flgtimecheck == 1 )  call system_clock(time_st, count_p_sec, count_max)
!
!
!
! Test two building types
!
  read(*,*)nbtype
  allocate(cases(nbtype))
  allocate(blds(nbtype))
!
!
!
  resolution = 1
!
!
!
  write(*,*)'++++ STEBBS ++++'
  write(*,*)'    + Input namelist case name...'
!
!
!
  do i = 1, nbtype, 1
      read(*,*)cases(i)
      write(*,*)'    + Namelist : ' , './BuildClasses/'// trim(cases(i))//'.nml'
!
      call create_building(cases(i),blds(i),i)
!
  enddo
!
!
!
  call readsuewsout()
!
!
!
  do i = 1, nbtype, 1
      call suewsstebbscouple(blds(i))
  enddo
!
!
!
  do i = 1, nbtype, 1
  do j = 1, 4, 1
    close(j + 100*i)
  enddo
  enddo
!
!
!
  if ( flgtimecheck == 1 ) then
      call system_clock(time_ed)
      write(*,*) '    + Elapsed time : ' , real( time_ed - time_st ) / real(count_p_sec)
  endif
!
!
!
#ifndef STEBBSONLINE
  end program
#else
  end subroutine
#endif
