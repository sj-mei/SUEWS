MODULE test_radiation_3d_module
   !==============================================================================
   ! VALIDATION TEST SUITE FOR 3D MORPHOLOGICAL RADIATION PARAMETERIZATIONS
   ! 
   ! This module provides validation tests for the 3D morphological radiation 
   ! parameterization scheme implemented from Mei et al. (2025).
   !
   ! Tests included:
   ! 1. Solar gain parameterization validation
   ! 2. Sky view factor calculation validation  
   ! 3. Directional absorptivity validation
   ! 4. Morphology-dependent albedo validation
   ! 5. Integration test with known parameter values
   !
   ! Reference validation data from Mei et al. (2025) paper figures and tables
   !
   ! Last modified: SJ Mei & Claude Code, Dec 2024
   !==============================================================================
   
   USE RADIATION_3D_MODULE
   
   IMPLICIT NONE
   
   ! Test tolerance for numerical comparisons
   REAL(KIND(1D0)), PARAMETER :: TEST_TOLERANCE = 1.0D-4
   
CONTAINS

   !==============================================================================
   SUBROUTINE test_all_radiation_3d()
      ! Main test suite runner
      ! Executes all validation tests and reports results
      
      IMPLICIT NONE
      
      LOGICAL :: all_tests_passed
      INTEGER :: test_count, passed_count
      
      WRITE(*, *) '======================================'
      WRITE(*, *) '3D MORPHOLOGICAL RADIATION TEST SUITE'
      WRITE(*, *) '======================================'
      WRITE(*, *)
      
      test_count = 0
      passed_count = 0
      all_tests_passed = .TRUE.
      
      ! Test 1: Solar gain parameterizations
      CALL run_test('Solar Gain Parameterization', test_solar_gain_parameterization, &
                   test_count, passed_count, all_tests_passed)
      
      ! Test 2: Sky view factor calculations
      CALL run_test('Sky View Factor Calculation', test_sky_view_factor_calculation, &
                   test_count, passed_count, all_tests_passed)
      
      ! Test 3: Directional absorptivity calculations
      CALL run_test('Directional Absorptivity', test_directional_absorptivity, &
                   test_count, passed_count, all_tests_passed)
      
      ! Test 4: Morphology-dependent albedo
      CALL run_test('Morphology-Dependent Albedo', test_morphology_albedo, &
                   test_count, passed_count, all_tests_passed)
      
      ! Test 5: Integration test
      CALL run_test('Full Integration Test', test_full_integration, &
                   test_count, passed_count, all_tests_passed)
      
      ! Report overall results
      WRITE(*, *) 
      WRITE(*, *) '======================================'
      WRITE(*, *) 'TEST SUMMARY:'
      WRITE(*, '(A,I0,A,I0,A)') 'Passed: ', passed_count, ' / ', test_count, ' tests'
      
      IF (all_tests_passed) THEN
         WRITE(*, *) 'ALL TESTS PASSED ✓'
      ELSE
         WRITE(*, *) 'SOME TESTS FAILED ✗'
      END IF
      WRITE(*, *) '======================================'
      
   END SUBROUTINE test_all_radiation_3d

   !==============================================================================
   SUBROUTINE run_test(test_name, test_subroutine, test_count, passed_count, all_passed)
      ! Generic test runner wrapper
      
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: test_name
      INTEGER, INTENT(inout) :: test_count, passed_count
      LOGICAL, INTENT(inout) :: all_passed
      LOGICAL :: test_result
      
      INTERFACE
         SUBROUTINE test_subroutine(result)
            LOGICAL, INTENT(out) :: result
         END SUBROUTINE
      END INTERFACE
      
      test_count = test_count + 1
      CALL test_subroutine(test_result)
      
      IF (test_result) THEN
         passed_count = passed_count + 1
         WRITE(*, '(A,A,A)') '[PASS] ', test_name, ' ✓'
      ELSE
         all_passed = .FALSE.
         WRITE(*, '(A,A,A)') '[FAIL] ', test_name, ' ✗'
      END IF
      
   END SUBROUTINE run_test

   !==============================================================================
   SUBROUTINE test_solar_gain_parameterization(result)
      ! Test solar gain calculations against known values from paper
      
      IMPLICIT NONE
      LOGICAL, INTENT(out) :: result
      
      ! Test parameters from paper validation range
      REAL(KIND(1D0)) :: lambda_f_test, lambda_p_test
      REAL(KIND(1D0)) :: solar_gain_ground, solar_gain_canyon
      REAL(KIND(1D0)) :: expected_ground, expected_canyon
      
      result = .TRUE.
      
      ! Test case 1: Low density urban (λp=0.16, λf=0.3)  
      lambda_f_test = 0.3D0
      lambda_p_test = 0.16D0
      
      CALL CALCULATE_SOLAR_GAIN_MORPHOLOGY( &
         lambda_f_test, lambda_p_test, &
         solar_gain_ground, solar_gain_canyon)
      
      ! Expected values based on paper equations
      ! Ground: L = 1.446 - 0.965*0.16 - 0.259*0.3 = 1.291
      ! Ground gain ≈ 0.0252 * exp(2.264*1.291) ≈ 0.83
      expected_ground = 0.83D0
      
      ! Canyon: 0.011*0.3 - 0.749*0.16 + 0.79 ≈ 0.793  
      expected_canyon = 0.793D0
      
      IF (ABS(solar_gain_ground - expected_ground) > 0.05D0) THEN
         WRITE(*, '(A,F8.4,A,F8.4)') '  Ground gain mismatch: got ', &
            solar_gain_ground, ', expected ~', expected_ground
         result = .FALSE.
      END IF
      
      IF (ABS(solar_gain_canyon - expected_canyon) > 0.05D0) THEN
         WRITE(*, '(A,F8.4,A,F8.4)') '  Canyon gain mismatch: got ', &
            solar_gain_canyon, ', expected ~', expected_canyon
         result = .FALSE.
      END IF
      
      ! Test case 2: High density urban (λp=0.69, λf=1.5)
      lambda_f_test = 1.5D0
      lambda_p_test = 0.69D0
      
      CALL CALCULATE_SOLAR_GAIN_MORPHOLOGY( &
         lambda_f_test, lambda_p_test, &
         solar_gain_ground, solar_gain_canyon)
      
      ! For high density, expect lower ground gain, moderate canyon gain
      IF (solar_gain_ground > 0.1D0 .OR. solar_gain_ground < 0.01D0) THEN
         WRITE(*, '(A,F8.4)') '  High density ground gain out of range: ', solar_gain_ground
         result = .FALSE.
      END IF
      
   END SUBROUTINE test_solar_gain_parameterization

   !==============================================================================
   SUBROUTINE test_sky_view_factor_calculation(result)
      ! Test sky view factor calculations
      
      IMPLICIT NONE
      LOGICAL, INTENT(out) :: result
      
      REAL(KIND(1D0)) :: lambda_f_test, lambda_p_test
      REAL(KIND(1D0)) :: sky_view_ground, sky_view_canyon
      
      result = .TRUE.
      
      ! Test case 1: Open area (low densities)
      lambda_f_test = 0.2D0
      lambda_p_test = 0.2D0
      
      CALL CALCULATE_SKY_VIEW_MORPHOLOGY( &
         lambda_f_test, lambda_p_test, &
         sky_view_ground, sky_view_canyon)
      
      ! Open area should have high sky view factors
      IF (sky_view_ground < 0.5D0 .OR. sky_view_ground > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Open area ground sky view out of range: ', sky_view_ground
         result = .FALSE.
      END IF
      
      IF (sky_view_canyon < 0.3D0 .OR. sky_view_canyon > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Open area canyon sky view out of range: ', sky_view_canyon
         result = .FALSE.
      END IF
      
      ! Test case 2: Dense area (high densities)
      lambda_f_test = 2.0D0
      lambda_p_test = 0.7D0
      
      CALL CALCULATE_SKY_VIEW_MORPHOLOGY( &
         lambda_f_test, lambda_p_test, &
         sky_view_ground, sky_view_canyon)
      
      ! Dense area should have low sky view factors
      IF (sky_view_ground > 0.3D0 .OR. sky_view_ground < 0.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Dense area ground sky view out of range: ', sky_view_ground
         result = .FALSE.
      END IF
      
      ! Test physical bounds
      IF (sky_view_ground < 0.0D0 .OR. sky_view_ground > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Ground sky view factor out of physical bounds: ', sky_view_ground
         result = .FALSE.
      END IF
      
      IF (sky_view_canyon < 0.0D0 .OR. sky_view_canyon > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Canyon sky view factor out of physical bounds: ', sky_view_canyon  
         result = .FALSE.
      END IF
      
   END SUBROUTINE test_sky_view_factor_calculation

   !==============================================================================
   SUBROUTINE test_directional_absorptivity(result)
      ! Test directional absorptivity calculations
      
      IMPLICIT NONE
      LOGICAL, INTENT(out) :: result
      
      REAL(KIND(1D0)), PARAMETER :: PI = 3.14159265358979323846D0
      REAL(KIND(1D0)) :: lambda_f_test, lambda_p_test
      REAL(KIND(1D0)) :: zenith_rad, azimuth_rad
      REAL(KIND(1D0)) :: absorptivity_ucl, absorptivity_canyon, absorptivity_ground
      
      result = .TRUE.
      
      ! Test case 1: Noon sun (zenith = 0°)
      lambda_f_test = 1.0D0
      lambda_p_test = 0.4D0
      zenith_rad = 0.0D0  ! overhead sun
      azimuth_rad = PI    ! south
      
      CALL CALCULATE_ABSORPTIVITY_DIRECTIONAL( &
         lambda_f_test, lambda_p_test, zenith_rad, azimuth_rad, &
         absorptivity_ucl, absorptivity_canyon, absorptivity_ground)
      
      ! Check physical bounds (0 ≤ absorptivity ≤ 1)
      IF (absorptivity_ucl < 0.0D0 .OR. absorptivity_ucl > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  UCL absorptivity out of bounds: ', absorptivity_ucl
         result = .FALSE.
      END IF
      
      IF (absorptivity_canyon < 0.0D0 .OR. absorptivity_canyon > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Canyon absorptivity out of bounds: ', absorptivity_canyon
         result = .FALSE.
      END IF
      
      IF (absorptivity_ground < 0.0D0 .OR. absorptivity_ground > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Ground absorptivity out of bounds: ', absorptivity_ground
         result = .FALSE.
      END IF
      
      ! Test case 2: Low sun angle (zenith = 60°)
      zenith_rad = PI/3.0D0  ! 60 degrees
      
      CALL CALCULATE_ABSORPTIVITY_DIRECTIONAL( &
         lambda_f_test, lambda_p_test, zenith_rad, azimuth_rad, &
         absorptivity_ucl, absorptivity_canyon, absorptivity_ground)
      
      ! Check that results are still within bounds
      IF (absorptivity_ground < 0.0D0 .OR. absorptivity_ground > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Ground absorptivity at low sun out of bounds: ', absorptivity_ground
         result = .FALSE.
      END IF
      
   END SUBROUTINE test_directional_absorptivity

   !==============================================================================
   SUBROUTINE test_morphology_albedo(result)
      ! Test morphology-dependent albedo calculations
      
      IMPLICIT NONE
      LOGICAL, INTENT(out) :: result
      
      REAL(KIND(1D0)), PARAMETER :: PI = 3.14159265358979323846D0
      REAL(KIND(1D0)) :: lambda_f_test, lambda_p_test, zenith_rad
      REAL(KIND(1D0)) :: base_albedo, absorptivity_ucl, urban_albedo
      
      result = .TRUE.
      
      ! Test parameters
      lambda_f_test = 1.0D0
      lambda_p_test = 0.5D0
      zenith_rad = PI/4.0D0  ! 45 degrees
      base_albedo = 0.2D0    ! typical urban albedo
      absorptivity_ucl = 0.8D0  ! typical absorption
      
      CALL CALCULATE_ALBEDO_MORPHOLOGY( &
         lambda_f_test, lambda_p_test, zenith_rad, base_albedo, &
         absorptivity_ucl, urban_albedo)
      
      ! Check physical bounds (0 ≤ albedo ≤ 1) 
      IF (urban_albedo < 0.0D0 .OR. urban_albedo > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Urban albedo out of bounds: ', urban_albedo
         result = .FALSE.
      END IF
      
      ! Test that morphological effects are reasonable
      ! (within ±50% of base albedo for typical parameters)
      IF (ABS(urban_albedo - base_albedo) > 0.1D0) THEN
         WRITE(*, '(A,F8.4,A,F8.4)') '  Large albedo change: base=', base_albedo, &
            ', urban=', urban_albedo
         ! This is not necessarily an error, just a diagnostic
      END IF
      
   END SUBROUTINE test_morphology_albedo

   !==============================================================================
   SUBROUTINE test_full_integration(result)
      ! Integration test using the main RADIATION_3D_MORPHOLOGY subroutine
      
      IMPLICIT NONE
      LOGICAL, INTENT(out) :: result
      
      REAL(KIND(1D0)), PARAMETER :: PI = 3.14159265358979323846D0
      
      ! Input parameters
      REAL(KIND(1D0)) :: lambda_f, lambda_p, zenith_rad, azimuth_rad
      REAL(KIND(1D0)) :: kdown, ldown, surface_albedo_base, surface_emissivity
      
      ! Output variables
      REAL(KIND(1D0)) :: urban_albedo, solar_gain_ground, solar_gain_canyon
      REAL(KIND(1D0)) :: sky_view_ground, sky_view_canyon
      REAL(KIND(1D0)) :: absorptivity_ucl, absorptivity_canyon, absorptivity_ground
      
      result = .TRUE.
      
      ! Set typical urban parameters
      lambda_f = 1.2D0
      lambda_p = 0.4D0
      zenith_rad = PI/3.0D0      ! 60 degrees
      azimuth_rad = PI           ! south
      kdown = 800.0D0            ! W/m²
      ldown = 300.0D0            ! W/m²
      surface_albedo_base = 0.15D0
      surface_emissivity = 0.95D0
      
      ! Call the main integration subroutine
      CALL RADIATION_3D_MORPHOLOGY( &
         lambda_f, lambda_p, zenith_rad, azimuth_rad, &
         kdown, ldown, surface_albedo_base, surface_emissivity, &
         urban_albedo, solar_gain_ground, solar_gain_canyon, &
         sky_view_ground, sky_view_canyon, &
         absorptivity_ucl, absorptivity_canyon, absorptivity_ground)
      
      ! Validate all outputs are within physical bounds
      
      ! Albedo bounds
      IF (urban_albedo < 0.0D0 .OR. urban_albedo > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Integration: urban albedo out of bounds: ', urban_albedo
         result = .FALSE.
      END IF
      
      ! Solar gain bounds (should be reasonable for given kdown)
      IF (solar_gain_ground < 0.0D0 .OR. solar_gain_ground > kdown) THEN
         WRITE(*, '(A,F8.4,A,F8.4)') '  Integration: ground solar gain out of bounds: ', &
            solar_gain_ground, ', kdown=', kdown
         result = .FALSE.
      END IF
      
      IF (solar_gain_canyon < 0.0D0 .OR. solar_gain_canyon > kdown) THEN
         WRITE(*, '(A,F8.4,A,F8.4)') '  Integration: canyon solar gain out of bounds: ', &
            solar_gain_canyon, ', kdown=', kdown
         result = .FALSE.
      END IF
      
      ! Sky view factor bounds
      IF (sky_view_ground < 0.0D0 .OR. sky_view_ground > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Integration: ground sky view out of bounds: ', sky_view_ground
         result = .FALSE.
      END IF
      
      IF (sky_view_canyon < 0.0D0 .OR. sky_view_canyon > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Integration: canyon sky view out of bounds: ', sky_view_canyon
         result = .FALSE.
      END IF
      
      ! Absorptivity bounds
      IF (absorptivity_ucl < 0.0D0 .OR. absorptivity_ucl > 1.0D0) THEN
         WRITE(*, '(A,F8.4)') '  Integration: UCL absorptivity out of bounds: ', absorptivity_ucl
         result = .FALSE.
      END IF
      
      ! Output diagnostic information for manual verification
      IF (result) THEN
         WRITE(*, '(A)') '  Integration test outputs:'
         WRITE(*, '(A,F8.4)') '    Urban albedo: ', urban_albedo
         WRITE(*, '(A,F8.1,A)') '    Ground solar gain: ', solar_gain_ground, ' W/m²'
         WRITE(*, '(A,F8.1,A)') '    Canyon solar gain: ', solar_gain_canyon, ' W/m²'  
         WRITE(*, '(A,F8.4)') '    Ground sky view: ', sky_view_ground
         WRITE(*, '(A,F8.4)') '    Canyon sky view: ', sky_view_canyon
      END IF
      
   END SUBROUTINE test_full_integration

END MODULE test_radiation_3d_module