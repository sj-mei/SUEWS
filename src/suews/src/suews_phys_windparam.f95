MODULE windparam_module
   ! CFD-based wind parameterization scheme
   ! Based on: Niu, J., Mei, S.-J., & Sun, T. (2025). Efficient city-scale wind mapping from building morphology: 
   !           A CFD-based parameterization scheme. Sustainable Cities and Society, 131, 106688.
   !
   ! This module implements the nonlinear regression-based wind speed parameterization
   ! that uses both plan area density (λp) and frontal area density (λf) to predict
   ! wind speeds at canopy and pedestrian levels in urban environments.
   !
   ! Last modified:
   ! SJ Mei & Claude Code, Dec 2024 - Initial implementation for SUEWS integration
   
   USE PhysConstants, ONLY: eps_fp
   IMPLICIT NONE

   ! Constants from the CFD parameterization paper
   REAL(KIND(1D0)), PARAMETER :: &
      ! Canopy layer L parameters: L = -0.36*λp + 0.98*λf + 0.1
      a_canopy = -0.36D0, &    ! λp coefficient for canopy layer
      b_canopy = 0.98D0, &     ! λf coefficient for canopy layer  
      c_canopy = 0.1D0, &      ! constant term for canopy layer
      
      ! Pedestrian level L parameters: L = -0.26*λp + 0.97*λf + 0.1
      a_pedestrian = -0.26D0, &  ! λp coefficient for pedestrian level
      b_pedestrian = 0.97D0, &   ! λf coefficient for pedestrian level
      c_pedestrian = 0.1D0, &    ! constant term for pedestrian level
      
      ! Wind speed equation parameters
      coeff_1_canopy = 0.017D0, &    ! U*c = 0.017*L^(-2.2) + 0.09
      coeff_2_canopy = -2.2D0, &     ! power coefficient for canopy
      coeff_3_canopy = 0.09D0, &     ! additive constant for canopy
      
      coeff_1_pedestrian = -43D0, &  ! U*p = exp(-43*L^4) * (0.74 - 2.52*L) / L + 0.297
      coeff_2_pedestrian = 4D0, &    ! power in exponential term
      coeff_3_pedestrian = 0.74D0, & ! coefficient in numerator
      coeff_4_pedestrian = -2.52D0, & ! coefficient for L term
      coeff_5_pedestrian = 0.297D0, & ! additive constant
      
      ! Morphological parameter bounds (from paper validation range)
      lambda_p_min = 0.16D0, &      ! minimum plan area density
      lambda_p_max = 0.69D0, &      ! maximum plan area density  
      lambda_f_min = 0.08D0, &      ! minimum frontal area density
      lambda_f_max = 1.74D0, &      ! maximum frontal area density (at 25m height)
      
      ! Physical constants
      kappa = 0.41D0                 ! von Karman constant

CONTAINS

   SUBROUTINE CFD_WindSpeed( &
      lambda_p, lambda_f, &           ! input: morphological parameters
      U_ref, z_ref, &                 ! input: reference wind speed and height
      building_height, &              ! input: average building height
      U_canopy, U_pedestrian)         ! output: wind speeds
      
      ! Calculates wind speeds at canopy and pedestrian levels using CFD-based parameterization
      ! 
      ! Input:
      !   lambda_p      - Plan area density (building footprint fraction) [-]
      !   lambda_f      - Frontal area density (frontal area per unit volume) [-]  
      !   U_ref         - Reference wind speed [m/s]
      !   z_ref         - Reference height [m]
      !   building_height - Average building height [m]
      ! Output:
      !   U_canopy      - Average wind speed in urban canopy layer [m/s]
      !   U_pedestrian  - Average wind speed at pedestrian level (2m) [m/s]
      
      IMPLICIT NONE
      
      ! Input variables
      REAL(KIND(1D0)), INTENT(in) :: lambda_p        ! plan area density [-]
      REAL(KIND(1D0)), INTENT(in) :: lambda_f        ! frontal area density [-]
      REAL(KIND(1D0)), INTENT(in) :: U_ref          ! reference wind speed [m/s]
      REAL(KIND(1D0)), INTENT(in) :: z_ref          ! reference height [m]
      REAL(KIND(1D0)), INTENT(in) :: building_height ! average building height [m]
      
      ! Output variables
      REAL(KIND(1D0)), INTENT(out) :: U_canopy      ! canopy layer wind speed [m/s]
      REAL(KIND(1D0)), INTENT(out) :: U_pedestrian  ! pedestrian level wind speed [m/s]
      
      ! Local variables
      REAL(KIND(1D0)) :: L_canopy           ! L parameter for canopy layer
      REAL(KIND(1D0)) :: L_pedestrian      ! L parameter for pedestrian level
      REAL(KIND(1D0)) :: U_star_canopy     ! normalized canopy wind speed
      REAL(KIND(1D0)) :: U_star_pedestrian ! normalized pedestrian wind speed
      REAL(KIND(1D0)) :: u_friction        ! friction velocity
      REAL(KIND(1D0)) :: reference_profile_integral ! integral of reference wind profile
      
      ! Validate input parameters
      CALL validate_morphological_parameters(lambda_p, lambda_f)
      
      ! Calculate L parameters for both levels
      L_canopy = a_canopy * lambda_p + b_canopy * lambda_f + c_canopy
      L_pedestrian = a_pedestrian * lambda_p + b_pedestrian * lambda_f + c_pedestrian
      
      ! Calculate normalized wind speeds using CFD parameterization equations
      CALL calculate_normalized_canopy_wind(L_canopy, U_star_canopy)
      CALL calculate_normalized_pedestrian_wind(lambda_f, lambda_p, U_star_pedestrian)
      
      ! Convert normalized wind speeds to physical wind speeds
      ! Based on paper normalization: U*c = AUc/∫U(z)dz and U*p = AUp/(Hu*)
      
      ! Calculate friction velocity from reference conditions
      u_friction = calculate_friction_velocity(U_ref, z_ref, building_height)
      
      ! Calculate reference wind profile integral (∫U(z)dz from 0 to H)
      reference_profile_integral = calculate_wind_profile_integral(u_friction, building_height)
      
      ! Denormalize wind speeds
      ! For canopy: U*c = A*Uc/∫U(z)dz  => Uc = U*c * ∫U(z)dz / A
      ! Assuming A ≈ building spacing, use building_height as characteristic length
      U_canopy = U_star_canopy * reference_profile_integral / building_height
      
      ! For pedestrian: Up = Up_star * U_star using correct formula
      ! where Up_star = (lambda_f / lambda_p) * ln(11/1) and U_star is friction velocity  
      U_pedestrian = U_star_pedestrian * u_friction
      
      ! Apply reasonable bounds to prevent unrealistic values
      U_canopy = MAX(0.1D0, MIN(U_canopy, 2.0D0 * U_ref))
      U_pedestrian = MAX(0.05D0, MIN(U_pedestrian, 1.5D0 * U_ref))
      
   END SUBROUTINE CFD_WindSpeed

   SUBROUTINE calculate_normalized_canopy_wind(L, U_star)
      ! Calculate normalized canopy wind speed: U*c = 0.017*L^(-2.2) + 0.09
      IMPLICIT NONE
      REAL(KIND(1D0)), INTENT(in) :: L
      REAL(KIND(1D0)), INTENT(out) :: U_star
      
      ! Handle the L^(-2.2) term safely
      IF (ABS(L) > eps_fp) THEN
         U_star = coeff_1_canopy * (L**coeff_2_canopy) + coeff_3_canopy
      ELSE
         ! For L ≈ 0, use limiting behavior
         U_star = coeff_3_canopy + 1.0D0  ! Large value when L approaches zero
      END IF
      
      ! Ensure positive values
      U_star = MAX(0.01D0, U_star)
      
   END SUBROUTINE calculate_normalized_canopy_wind

   SUBROUTINE calculate_normalized_pedestrian_wind(lambda_f, lambda_p, Up_star)
      ! Calculate Up_star using: Up_star = (lambda_f / lambda_p) * ln(11/1)
      ! Based on CFD parameterization approach
      IMPLICIT NONE
      REAL(KIND(1D0)), INTENT(in) :: lambda_f  ! frontal area density [-]
      REAL(KIND(1D0)), INTENT(in) :: lambda_p  ! plan area density [-]
      REAL(KIND(1D0)), INTENT(out) :: Up_star  ! normalized pedestrian wind coefficient
      
      REAL(KIND(1D0)) :: ln_ratio
      
      ! Calculate ln(11/1) = ln(11) ≈ 2.398
      ln_ratio = LOG(11.0D0)
      
      ! Calculate Up_star = (lambda_f / lambda_p) * ln(11/1)
      IF (ABS(lambda_p) > eps_fp) THEN
         Up_star = (lambda_f / lambda_p) * ln_ratio
      ELSE
         ! If lambda_p is nearly zero, use a default value
         Up_star = lambda_f * ln_ratio  ! Assume lambda_p = 1 as fallback
      END IF
      
      ! Ensure reasonable values
      Up_star = MAX(0.01D0, MIN(Up_star, 10.0D0))
      
   END SUBROUTINE calculate_normalized_pedestrian_wind

   FUNCTION calculate_friction_velocity(U_ref, z_ref, z_canopy) RESULT(u_star)
      ! Calculate friction velocity from reference wind conditions
      ! Using logarithmic wind profile assumption
      IMPLICIT NONE
      REAL(KIND(1D0)), INTENT(in) :: U_ref, z_ref, z_canopy
      REAL(KIND(1D0)) :: u_star
      
      REAL(KIND(1D0)) :: z0, zd
      
      ! Estimate roughness parameters based on canopy height
      ! Using rule of thumb: z0 ≈ 0.1*H, zd ≈ 0.7*H for urban areas
      z0 = 0.1D0 * z_canopy
      zd = 0.7D0 * z_canopy
      
      ! Calculate friction velocity: u* = κ*U_ref / ln((z_ref-zd)/z0)
      IF (z_ref > zd + z0) THEN
         u_star = kappa * U_ref / LOG((z_ref - zd) / z0)
      ELSE
         ! Fallback for problematic reference heights
         u_star = 0.1D0 * U_ref  ! Simple approximation
      END IF
      
      ! Ensure reasonable values
      u_star = MAX(0.01D0, MIN(u_star, 2.0D0))
      
   END FUNCTION calculate_friction_velocity

   FUNCTION calculate_wind_profile_integral(u_star, H) RESULT(integral)
      ! Calculate integral of wind profile from ground to canopy height
      ! ∫[0 to H] U(z) dz using logarithmic profile
      IMPLICIT NONE
      REAL(KIND(1D0)), INTENT(in) :: u_star, H
      REAL(KIND(1D0)) :: integral
      
      REAL(KIND(1D0)) :: z0, zd
      
      ! Roughness parameters
      z0 = 0.1D0 * H
      zd = 0.7D0 * H
      
      ! Analytical integration of logarithmic profile
      ! This is a simplified approximation - exact integration is more complex
      integral = (u_star / kappa) * H * LOG(H / z0)
      
      ! Ensure positive value
      integral = MAX(0.1D0 * H, integral)
      
   END FUNCTION calculate_wind_profile_integral

   SUBROUTINE validate_morphological_parameters(lambda_p, lambda_f)
      ! Validate that morphological parameters are within acceptable ranges
      IMPLICIT NONE
      REAL(KIND(1D0)), INTENT(in) :: lambda_p, lambda_f
      
      IF (lambda_p < 0.0D0 .OR. lambda_p > 1.0D0) THEN
         CALL errorHint(45, 'CFD wind param: Plan area density out of range [0,1]', lambda_p, -999.0D0, -999)
      END IF
      
      IF (lambda_f < 0.0D0 .OR. lambda_f > 2.0D0) THEN
         CALL errorHint(45, 'CFD wind param: Frontal area density out of range [0,2]', lambda_f, -999.0D0, -999)
      END IF
      
      IF (lambda_p < lambda_p_min .OR. lambda_p > lambda_p_max) THEN
         CALL errorHint(46, 'CFD wind param: λp outside validated range [0.16,0.69], results may be unreliable', &
                        lambda_p, -999.0D0, -999)
      END IF
      
      IF (lambda_f < lambda_f_min .OR. lambda_f > lambda_f_max) THEN
         CALL errorHint(46, 'CFD wind param: λf outside validated range [0.08,1.74], results may be unreliable', &
                        lambda_f, -999.0D0, -999)
      END IF
      
   END SUBROUTINE validate_morphological_parameters

   SUBROUTINE CFD_WindProfile( &
      DiagMethod, &
      lambda_p, lambda_f, &           ! input: morphological parameters
      U_ref, z_ref, &                 ! input: reference conditions
      building_height, &              ! input: building characteristics
      z_levels, n_levels, &           ! input: height levels for profile calculation
      U_profile)                      ! output: wind speed profile
      
      ! Calculate wind speed profile at specified height levels
      ! Integrates with SUEWS diagnostic framework
      !
      ! This subroutine provides multi-level wind speed calculations for integration
      ! with existing SUEWS RSL profile capabilities and diagnostic output
      
      IMPLICIT NONE
      
      ! Input variables
      INTEGER, INTENT(in) :: DiagMethod                      ! diagnostic method selection
      REAL(KIND(1D0)), INTENT(in) :: lambda_p, lambda_f     ! morphological parameters
      REAL(KIND(1D0)), INTENT(in) :: U_ref, z_ref           ! reference conditions
      REAL(KIND(1D0)), INTENT(in) :: building_height        ! average building height
      INTEGER, INTENT(in) :: n_levels                       ! number of height levels
      REAL(KIND(1D0)), DIMENSION(n_levels), INTENT(in) :: z_levels  ! height levels [m]
      
      ! Output variables  
      REAL(KIND(1D0)), DIMENSION(n_levels), INTENT(out) :: U_profile ! wind speeds [m/s]
      
      ! Local variables
      REAL(KIND(1D0)) :: U_canopy, U_pedestrian            ! base wind speeds
      REAL(KIND(1D0)) :: u_friction                        ! friction velocity
      INTEGER :: i
      
      ! Calculate base wind speeds using CFD parameterization
      CALL CFD_WindSpeed(lambda_p, lambda_f, U_ref, z_ref, building_height, &
                        U_canopy, U_pedestrian)
      
      ! Calculate friction velocity for profile extrapolation
      u_friction = calculate_friction_velocity(U_ref, z_ref, building_height)
      
      ! Calculate wind speed at each level
      DO i = 1, n_levels
         IF (z_levels(i) <= 2.0D0) THEN
            ! Use pedestrian-level parameterization for heights ≤ 2m
            U_profile(i) = U_pedestrian * (z_levels(i) / 2.0D0)**0.1D0
         ELSE IF (z_levels(i) <= building_height) THEN
            ! Use canopy parameterization within urban canopy layer
            U_profile(i) = U_canopy * (z_levels(i) / building_height)**0.3D0
         ELSE
            ! Above canopy: use logarithmic profile from canopy-top wind speed
            U_profile(i) = extrapolate_above_canopy(U_canopy, building_height, z_levels(i), u_friction)
         END IF
         
         ! Ensure reasonable bounds
         U_profile(i) = MAX(0.1D0, MIN(U_profile(i), 3.0D0 * U_ref))
      END DO
      
   END SUBROUTINE CFD_WindProfile

   FUNCTION extrapolate_above_canopy(U_canopy, H, z, u_star) RESULT(U_z)
      ! Extrapolate wind speed above canopy using logarithmic profile
      IMPLICIT NONE
      REAL(KIND(1D0)), INTENT(in) :: U_canopy, H, z, u_star
      REAL(KIND(1D0)) :: U_z
      
      REAL(KIND(1D0)) :: z0, zd
      
      ! Roughness parameters
      z0 = 0.1D0 * H
      zd = 0.7D0 * H
      
      ! Logarithmic profile above canopy
      IF (z > zd + z0 .AND. H > zd + z0) THEN
         U_z = U_canopy * LOG((z - zd) / z0) / LOG((H - zd) / z0)
      ELSE
         ! Simple power law as fallback
         U_z = U_canopy * (z / H)**0.2D0
      END IF
      
   END FUNCTION extrapolate_above_canopy

END MODULE windparam_module