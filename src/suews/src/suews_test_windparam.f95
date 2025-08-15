PROGRAM test_cfd_wind_parameterization
   ! Test program for CFD-based wind parameterization
   ! This program validates the implementation against values from the paper
   !
   ! Based on: Niu, J., Mei, S.-J., & Sun, T. (2025). Efficient city-scale wind mapping from building morphology: 
   !           A CFD-based parameterization scheme. Sustainable Cities and Society, 131, 106688.
   !
   ! Test cases validate:
   ! 1. L parameter calculation for different λp and λf values
   ! 2. Wind speed parameterization equations 
   ! 3. Morphological parameter conversions
   ! 4. Edge cases and boundary conditions
   !
   ! Last modified:
   ! SJ Mei & Claude Code, Dec 2024 - Initial validation test suite

   USE windparam_module
   USE morphconv_module
   IMPLICIT NONE

   ! Test parameters
   LOGICAL :: all_tests_passed = .TRUE.
   INTEGER :: test_count = 0
   INTEGER :: passed_count = 0

   WRITE(*,*) '=========================================='
   WRITE(*,*) 'CFD Wind Parameterization Validation Tests'
   WRITE(*,*) '=========================================='
   WRITE(*,*) ''

   ! Test 1: Basic parameter validation 
   CALL test_morphological_validation()
   
   ! Test 2: L parameter calculations
   CALL test_L_parameter_calculations()
   
   ! Test 3: Wind speed equations
   CALL test_wind_speed_equations()
   
   ! Test 4: Morphological conversions
   CALL test_morphological_conversions()
   
   ! Test 5: Integration test
   CALL test_cfd_wind_integration()

   ! Summary
   WRITE(*,*) ''
   WRITE(*,*) '=========================================='
   WRITE(*,'(A,I0,A,I0,A)') 'SUMMARY: ', passed_count, '/', test_count, ' tests passed'
   IF (all_tests_passed) THEN
      WRITE(*,*) 'All tests PASSED ✓'
      WRITE(*,*) 'CFD wind parameterization implementation is valid!'
   ELSE
      WRITE(*,*) 'Some tests FAILED ✗'
      WRITE(*,*) 'Please check the implementation.'
   END IF
   WRITE(*,*) '=========================================='

CONTAINS

   SUBROUTINE test_morphological_validation()
      WRITE(*,*) 'Test 1: Morphological Parameter Validation'
      WRITE(*,*) '-------------------------------------------'
      
      ! This test would normally call validation routines
      ! For now, we just check that the constants are defined correctly
      IF (lambda_p_min == 0.16D0 .AND. lambda_p_max == 0.69D0 .AND. &
          lambda_f_min == 0.08D0 .AND. lambda_f_max == 1.74D0) THEN
         WRITE(*,*) '  ✓ Parameter bounds correctly defined'
         CALL record_test_result(.TRUE.)
      ELSE
         WRITE(*,*) '  ✗ Parameter bounds incorrectly defined'  
         CALL record_test_result(.FALSE.)
      END IF
      WRITE(*,*) ''
   END SUBROUTINE test_morphological_validation

   SUBROUTINE test_L_parameter_calculations()
      REAL(KIND(1D0)) :: lambda_p, lambda_f, L_expected, L_calculated
      REAL(KIND(1D0)), PARAMETER :: tolerance = 1.0E-6
      
      WRITE(*,*) 'Test 2: L Parameter Calculations'
      WRITE(*,*) '--------------------------------'
      
      ! Test case from paper: typical urban morphology
      lambda_p = 0.25D0  ! 25% plan area density
      lambda_f = 0.5D0   ! 0.5 m⁻¹ frontal area density
      
      ! Expected L for canopy: L = -0.36*0.25 + 0.98*0.5 + 0.1 = 0.5
      L_expected = -0.36D0 * lambda_p + 0.98D0 * lambda_f + 0.1D0
      L_calculated = a_canopy * lambda_p + b_canopy * lambda_f + c_canopy
      
      IF (ABS(L_calculated - L_expected) < tolerance) THEN
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✓ Canopy L calculation: ', L_calculated, ' (expected: ', L_expected, ')'
         CALL record_test_result(.TRUE.)
      ELSE
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✗ Canopy L calculation: ', L_calculated, ' (expected: ', L_expected, ')'
         CALL record_test_result(.FALSE.)
      END IF
      
      ! Test pedestrian L parameter
      L_expected = -0.26D0 * lambda_p + 0.97D0 * lambda_f + 0.1D0
      L_calculated = a_pedestrian * lambda_p + b_pedestrian * lambda_f + c_pedestrian
      
      IF (ABS(L_calculated - L_expected) < tolerance) THEN
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✓ Pedestrian L calculation: ', L_calculated, ' (expected: ', L_expected, ')'
         CALL record_test_result(.TRUE.)
      ELSE
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✗ Pedestrian L calculation: ', L_calculated, ' (expected: ', L_expected, ')'
         CALL record_test_result(.FALSE.)
      END IF
      WRITE(*,*) ''
   END SUBROUTINE test_L_parameter_calculations

   SUBROUTINE test_wind_speed_equations()
      REAL(KIND(1D0)) :: L, U_star_canopy, U_star_pedestrian
      REAL(KIND(1D0)) :: expected_canopy, expected_pedestrian
      REAL(KIND(1D0)), PARAMETER :: tolerance = 1.0E-3
      
      WRITE(*,*) 'Test 3: Wind Speed Equations'
      WRITE(*,*) '----------------------------'
      
      ! Test case: L = 0.5 (typical urban value)
      L = 0.5D0
      
      ! Expected canopy wind speed: U*c = 0.017 * L^(-2.2) + 0.09
      expected_canopy = 0.017D0 * (L**(-2.2D0)) + 0.09D0
      CALL calculate_normalized_canopy_wind(L, U_star_canopy)
      
      IF (ABS(U_star_canopy - expected_canopy) < tolerance) THEN
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✓ Canopy wind equation: ', U_star_canopy, ' (expected: ', expected_canopy, ')'
         CALL record_test_result(.TRUE.)
      ELSE
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✗ Canopy wind equation: ', U_star_canopy, ' (expected: ', expected_canopy, ')'
         CALL record_test_result(.FALSE.)
      END IF
      
      ! Expected pedestrian wind speed: U*p = exp(-43*L^4) * (0.74 - 2.52*L) / L + 0.297
      expected_pedestrian = EXP(-43D0 * L**4) * (0.74D0 - 2.52D0 * L) / L + 0.297D0
      CALL calculate_normalized_pedestrian_wind(L, U_star_pedestrian)
      
      IF (ABS(U_star_pedestrian - expected_pedestrian) < tolerance) THEN
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✓ Pedestrian wind equation: ', U_star_pedestrian, ' (expected: ', expected_pedestrian, ')'
         CALL record_test_result(.TRUE.)
      ELSE
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✗ Pedestrian wind equation: ', U_star_pedestrian, ' (expected: ', expected_pedestrian, ')'
         CALL record_test_result(.FALSE.)
      END IF
      WRITE(*,*) ''
   END SUBROUTINE test_wind_speed_equations

   SUBROUTINE test_morphological_conversions()
      REAL(KIND(1D0)) :: FAI, PAI, building_height, surface_area
      REAL(KIND(1D0)) :: lambda_f, lambda_p
      REAL(KIND(1D0)) :: expected_lambda_f, expected_lambda_p
      REAL(KIND(1D0)), PARAMETER :: tolerance = 1.0E-3
      
      WRITE(*,*) 'Test 4: Morphological Parameter Conversions'
      WRITE(*,*) '-------------------------------------------'
      
      ! Test typical urban values
      FAI = 1.0D0           ! Frontal Area Index
      PAI = 0.25D0          ! Plan Area Index  
      building_height = 20.0D0  ! 20m average building height
      surface_area = 250000.0D0  ! 500m x 500m grid
      
      ! Expected conversions
      expected_lambda_p = PAI  ! λp ≈ PAI
      expected_lambda_f = FAI / building_height  ! λf ≈ FAI / height
      
      CALL convert_FAI_PAI_to_lambda(FAI, PAI, building_height, surface_area, lambda_f, lambda_p)
      
      IF (ABS(lambda_p - expected_lambda_p) < tolerance) THEN
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✓ λp conversion: ', lambda_p, ' (expected: ', expected_lambda_p, ')'
         CALL record_test_result(.TRUE.)
      ELSE
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✗ λp conversion: ', lambda_p, ' (expected: ', expected_lambda_p, ')'
         CALL record_test_result(.FALSE.)
      END IF
      
      IF (ABS(lambda_f - expected_lambda_f) < tolerance) THEN
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✓ λf conversion: ', lambda_f, ' (expected: ', expected_lambda_f, ')'
         CALL record_test_result(.TRUE.)
      ELSE
         WRITE(*,'(A,F6.3,A,F6.3)') '  ✗ λf conversion: ', lambda_f, ' (expected: ', expected_lambda_f, ')'
         CALL record_test_result(.FALSE.)
      END IF
      WRITE(*,*) ''
   END SUBROUTINE test_morphological_conversions

   SUBROUTINE test_cfd_wind_integration()
      REAL(KIND(1D0)) :: lambda_p, lambda_f, U_ref, z_ref, building_height
      REAL(KIND(1D0)) :: U_canopy, U_pedestrian
      
      WRITE(*,*) 'Test 5: CFD Wind Integration Test'
      WRITE(*,*) '--------------------------------'
      
      ! Test realistic urban scenario
      lambda_p = 0.3D0      ! 30% building coverage
      lambda_f = 0.6D0      ! 0.6 m⁻¹ frontal area density  
      U_ref = 5.0D0         ! 5 m/s reference wind speed
      z_ref = 50.0D0        ! 50m reference height
      building_height = 15.0D0  ! 15m average building height
      
      CALL CFD_WindSpeed(lambda_p, lambda_f, U_ref, z_ref, building_height, U_canopy, U_pedestrian)
      
      ! Basic sanity checks
      IF (U_canopy > 0.0D0 .AND. U_canopy < U_ref) THEN
         WRITE(*,'(A,F5.2,A)') '  ✓ Canopy wind speed: ', U_canopy, ' m/s (reasonable)'
         CALL record_test_result(.TRUE.)
      ELSE
         WRITE(*,'(A,F5.2,A)') '  ✗ Canopy wind speed: ', U_canopy, ' m/s (unreasonable)'
         CALL record_test_result(.FALSE.)
      END IF
      
      IF (U_pedestrian > 0.0D0 .AND. U_pedestrian < U_canopy) THEN
         WRITE(*,'(A,F5.2,A)') '  ✓ Pedestrian wind speed: ', U_pedestrian, ' m/s (reasonable)'
         CALL record_test_result(.TRUE.)
      ELSE
         WRITE(*,'(A,F5.2,A)') '  ✗ Pedestrian wind speed: ', U_pedestrian, ' m/s (unreasonable)'
         CALL record_test_result(.FALSE.)
      END IF
      
      ! Check that pedestrian wind is typically lower than canopy wind in urban areas
      IF (U_pedestrian < U_canopy) THEN
         WRITE(*,*) '  ✓ Pedestrian wind < canopy wind (expected for urban areas)'
         CALL record_test_result(.TRUE.)
      ELSE
         WRITE(*,*) '  ✗ Pedestrian wind >= canopy wind (unexpected for urban areas)'
         CALL record_test_result(.FALSE.)
      END IF
      WRITE(*,*) ''
   END SUBROUTINE test_cfd_wind_integration

   SUBROUTINE record_test_result(passed)
      LOGICAL, INTENT(in) :: passed
      test_count = test_count + 1
      IF (passed) THEN
         passed_count = passed_count + 1
      ELSE
         all_tests_passed = .FALSE.
      END IF
   END SUBROUTINE record_test_result

END PROGRAM test_cfd_wind_parameterization