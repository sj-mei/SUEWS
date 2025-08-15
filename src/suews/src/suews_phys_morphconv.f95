MODULE morphconv_module
   ! Morphological parameter conversion utilities
   ! Converts between SUEWS internal parameters (FAI, PAI) and 
   ! CFD wind parameterization parameters (λf, λp)
   !
   ! This module provides the interface between SUEWS existing morphological
   ! calculations and the CFD-based wind parameterization scheme requirements.
   !
   ! Last modified:
   ! SJ Mei & Claude Code, Dec 2024 - Initial implementation

   USE PhysConstants, ONLY: eps_fp
   IMPLICIT NONE

   ! Conversion constants and parameters
   REAL(KIND(1D0)), PARAMETER :: &
      default_building_width = 10.0D0, &      ! Default building width for area calculations [m]
      default_building_spacing = 20.0D0, &    ! Default building spacing [m]
      height_normalization = 10.0D0           ! Reference height for normalization [m]

CONTAINS

   SUBROUTINE convert_FAI_PAI_to_lambda( &
      FAI, PAI, &                           ! input: SUEWS morphological indices
      building_height, surface_area, &      ! input: geometric parameters
      lambda_f, lambda_p)                   ! output: CFD parameterization densities
      
      ! Converts SUEWS Frontal Area Index (FAI) and Plan Area Index (PAI) 
      ! to λf (frontal area density) and λp (plan area density) for CFD wind parameterization
      !
      ! Input:
      !   FAI              - Frontal Area Index from SUEWS [dimensionless]
      !   PAI              - Plan Area Index from SUEWS [dimensionless]  
      !   building_height  - Average building height [m]
      !   surface_area     - Total surface area [m²]
      ! Output:
      !   lambda_f         - Frontal area density (frontal area per unit volume) [m⁻¹]
      !   lambda_p         - Plan area density (building footprint fraction) [dimensionless]
      
      IMPLICIT NONE
      
      ! Input variables
      REAL(KIND(1D0)), INTENT(in) :: FAI              ! Frontal Area Index [-]
      REAL(KIND(1D0)), INTENT(in) :: PAI              ! Plan Area Index [-] 
      REAL(KIND(1D0)), INTENT(in) :: building_height  ! Average building height [m]
      REAL(KIND(1D0)), INTENT(in) :: surface_area     ! Surface area [m²]
      
      ! Output variables
      REAL(KIND(1D0)), INTENT(out) :: lambda_f        ! Frontal area density [m⁻¹]
      REAL(KIND(1D0)), INTENT(out) :: lambda_p        ! Plan area density [-]
      
      ! Local variables
      REAL(KIND(1D0)) :: volume_per_area               ! Volume per unit ground area [m]
      REAL(KIND(1D0)) :: effective_height              ! Effective height for calculations [m]
      
      ! Validate input parameters
      IF (FAI < 0.0D0 .OR. PAI < 0.0D0) THEN
         CALL errorHint(50, 'Morphological conversion: Negative FAI or PAI values detected', FAI, PAI, -999)
         lambda_f = 0.0D0
         lambda_p = 0.0D0
         RETURN
      END IF
      
      IF (building_height < eps_fp) THEN
         CALL errorHint(51, 'Morphological conversion: Building height too small, using default', &
                        building_height, height_normalization, -999)
         effective_height = height_normalization
      ELSE
         effective_height = building_height
      END IF
      
      ! Convert Plan Area Index to plan area density
      ! PAI in SUEWS represents the fraction of plan area occupied by buildings and trees
      ! λp is defined as the ratio of building plan area to total plan area
      ! For urban areas, this is approximately equivalent
      lambda_p = PAI
      
      ! Convert Frontal Area Index to frontal area density  
      ! FAI represents total frontal area of obstacles per unit ground area
      ! λf is frontal area per unit volume, so we need to divide by effective height
      ! λf = Frontal Area / (Ground Area × Height) = FAI / Height
      lambda_f = FAI / effective_height
      
      ! Apply bounds based on physical constraints and validation range from paper
      lambda_p = MIN(MAX(lambda_p, 0.0D0), 0.95D0)  ! Keep below 95% to avoid unrealistic density
      lambda_f = MIN(MAX(lambda_f, 0.0D0), 2.0D0)   ! Reasonable upper bound for urban areas
      
      ! Issue warnings if outside validated range from CFD paper
      IF (lambda_p < 0.16D0 .OR. lambda_p > 0.69D0) THEN
         CALL errorHint(52, 'λp outside CFD validation range [0.16,0.69], results may be less reliable', &
                        lambda_p, -999.0D0, -999)
      END IF
      
      IF (lambda_f < 0.08D0 .OR. lambda_f > 1.74D0) THEN
         CALL errorHint(52, 'λf outside CFD validation range [0.08,1.74], results may be less reliable', &
                        lambda_f, -999.0D0, -999)
      END IF
      
   END SUBROUTINE convert_FAI_PAI_to_lambda

   SUBROUTINE calculate_lambda_from_urban_geometry( &
      building_fraction, building_height, &        ! input: basic urban parameters
      average_building_width, average_spacing, &   ! input: geometric details
      lambda_f, lambda_p)                          ! output: CFD densities
      
      ! Direct calculation of λf and λp from basic urban geometric parameters
      ! This provides an alternative to FAI/PAI conversion for cases where
      ! detailed geometric information is available
      !
      ! Input:
      !   building_fraction     - Fraction of area covered by buildings [-]
      !   building_height       - Average building height [m]
      !   average_building_width - Average building width [m]
      !   average_spacing       - Average spacing between buildings [m]
      ! Output:
      !   lambda_f              - Frontal area density [m⁻¹]
      !   lambda_p              - Plan area density [-]
      
      IMPLICIT NONE
      
      ! Input variables
      REAL(KIND(1D0)), INTENT(in) :: building_fraction        ! Building area fraction [-]
      REAL(KIND(1D0)), INTENT(in) :: building_height          ! Average building height [m]
      REAL(KIND(1D0)), INTENT(in) :: average_building_width   ! Average building width [m]
      REAL(KIND(1D0)), INTENT(in) :: average_spacing          ! Average building spacing [m]
      
      ! Output variables
      REAL(KIND(1D0)), INTENT(out) :: lambda_f                ! Frontal area density [m⁻¹]
      REAL(KIND(1D0)), INTENT(out) :: lambda_p                ! Plan area density [-]
      
      ! Local variables
      REAL(KIND(1D0)) :: frontal_area_per_building             ! Frontal area of typical building [m²]
      REAL(KIND(1D0)) :: plan_area_per_building                ! Plan area of typical building [m²]
      REAL(KIND(1D0)) :: volume_per_building                   ! Volume occupied per building [m³]
      REAL(KIND(1D0)) :: total_area_per_building               ! Total area per building including spacing [m²]
      
      ! Validate inputs
      IF (building_fraction < 0.0D0 .OR. building_fraction > 1.0D0) THEN
         CALL errorHint(53, 'Invalid building fraction, must be [0,1]', building_fraction, -999.0D0, -999)
         lambda_f = 0.0D0
         lambda_p = 0.0D0
         RETURN
      END IF
      
      IF (building_height < eps_fp .OR. average_building_width < eps_fp .OR. average_spacing < eps_fp) THEN
         CALL errorHint(54, 'Invalid geometric parameters, using defaults', &
                        building_height, average_building_width, -999)
         ! Use default values
         CALL convert_building_fraction_to_lambda(building_fraction, building_height, lambda_f, lambda_p)
         RETURN
      END IF
      
      ! Calculate geometric properties
      frontal_area_per_building = building_height * average_building_width
      plan_area_per_building = average_building_width * average_building_width  ! Assume square buildings
      total_area_per_building = (average_building_width + average_spacing)**2   ! Include spacing
      volume_per_building = total_area_per_building * building_height
      
      ! Calculate densities
      ! λp = (Plan area of building) / (Total plan area) = building_fraction
      lambda_p = building_fraction
      
      ! λf = (Frontal area) / (Volume) = (Frontal area per building) / (Volume per building)
      lambda_f = frontal_area_per_building / volume_per_building
      
      ! Apply reasonable bounds
      lambda_p = MIN(MAX(lambda_p, 0.0D0), 0.95D0)
      lambda_f = MIN(MAX(lambda_f, 0.0D0), 2.0D0)
      
   END SUBROUTINE calculate_lambda_from_urban_geometry

   SUBROUTINE convert_building_fraction_to_lambda( &
      building_fraction, building_height, &         ! input: basic parameters
      lambda_f, lambda_p)                           ! output: CFD densities
      
      ! Simple conversion using building fraction and height with default assumptions
      ! This is a fallback method when detailed geometric information is not available
      !
      ! Input:
      !   building_fraction - Fraction of area covered by buildings [-]
      !   building_height   - Average building height [m]
      ! Output:
      !   lambda_f          - Frontal area density [m⁻¹] 
      !   lambda_p          - Plan area density [-]
      
      IMPLICIT NONE
      
      ! Input variables
      REAL(KIND(1D0)), INTENT(in) :: building_fraction    ! Building area fraction [-]
      REAL(KIND(1D0)), INTENT(in) :: building_height      ! Average building height [m]
      
      ! Output variables  
      REAL(KIND(1D0)), INTENT(out) :: lambda_f            ! Frontal area density [m⁻¹]
      REAL(KIND(1D0)), INTENT(out) :: lambda_p            ! Plan area density [-]
      
      ! Plan area density is directly the building fraction
      lambda_p = building_fraction
      
      ! For frontal area density, use empirical relationship based on typical urban geometries
      ! Assuming buildings have aspect ratio ~1 and are arranged in regular pattern
      IF (building_height > eps_fp) THEN
         ! λf ≈ √(λp) * H / (typical building dimension)
         ! This is a simplified approximation based on urban morphology studies
         lambda_f = SQRT(building_fraction) * building_height / default_building_spacing
      ELSE
         lambda_f = 0.0D0
      END IF
      
      ! Apply bounds
      lambda_p = MIN(MAX(lambda_p, 0.0D0), 0.95D0)
      lambda_f = MIN(MAX(lambda_f, 0.0D0), 2.0D0)
      
   END SUBROUTINE convert_building_fraction_to_lambda

   SUBROUTINE estimate_wind_exposure_adjustment( &
      lambda_p, lambda_f, wind_direction, &        ! input: morphological and wind parameters
      building_orientation_angle, &                ! input: building orientation relative to north
      lambda_f_adjusted)                           ! output: direction-adjusted frontal density
      
      ! Adjust frontal area density based on wind direction relative to building orientation
      ! This accounts for the directional nature of frontal area exposure
      !
      ! Input:
      !   lambda_p                  - Plan area density [-]
      !   lambda_f                  - Base frontal area density [m⁻¹] 
      !   wind_direction            - Wind direction [degrees from north]
      !   building_orientation_angle - Average building orientation [degrees from north]
      ! Output:
      !   lambda_f_adjusted         - Direction-adjusted frontal area density [m⁻¹]
      
      IMPLICIT NONE
      
      ! Input variables
      REAL(KIND(1D0)), INTENT(in) :: lambda_p                     ! Plan area density [-]
      REAL(KIND(1D0)), INTENT(in) :: lambda_f                     ! Base frontal area density [m⁻¹]
      REAL(KIND(1D0)), INTENT(in) :: wind_direction               ! Wind direction [degrees]
      REAL(KIND(1D0)), INTENT(in) :: building_orientation_angle   ! Building orientation [degrees]
      
      ! Output variables
      REAL(KIND(1D0)), INTENT(out) :: lambda_f_adjusted           ! Adjusted frontal area density [m⁻¹]
      
      ! Local variables
      REAL(KIND(1D0)) :: relative_angle      ! Wind direction relative to buildings [degrees]
      REAL(KIND(1D0)) :: exposure_factor     ! Factor accounting for directional exposure [-]
      REAL(KIND(1D0)), PARAMETER :: pi = 3.14159265358979323846D0
      
      ! Calculate relative wind direction
      relative_angle = ABS(wind_direction - building_orientation_angle)
      
      ! Normalize to [0, 90] degrees (considering symmetry)
      DO WHILE (relative_angle > 90.0D0)
         IF (relative_angle > 180.0D0) THEN
            relative_angle = 360.0D0 - relative_angle
         ELSE
            relative_angle = 180.0D0 - relative_angle
         END IF
      END DO
      
      ! Calculate exposure factor
      ! Maximum exposure when wind is perpendicular to buildings (90°)
      ! Minimum exposure when wind is parallel to buildings (0°)
      exposure_factor = SIN(relative_angle * pi / 180.0D0)
      
      ! Adjust frontal area density
      lambda_f_adjusted = lambda_f * exposure_factor
      
      ! Ensure non-negative result
      lambda_f_adjusted = MAX(0.0D0, lambda_f_adjusted)
      
   END SUBROUTINE estimate_wind_exposure_adjustment

END MODULE morphconv_module