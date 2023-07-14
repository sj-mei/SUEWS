MODULE SUEWS_HydroHeat_DTS
   USE allocateArray, ONLY: nsurf, nvegsurf

   PUBLIC :: HYDRO_STATE, HEAT_STATE
   PUBLIC :: allocate_hydro_state, dealloc_hydro_state, allocate_heat_state, dealloc_heat_state

   TYPE, PUBLIC :: HYDRO_STATE
      ! REAL(KIND(1D0)) :: runofftowater   ! Fraction of above-ground runoff flowing to water surface during flooding
      REAL(KIND(1D0)), DIMENSION(nsurf) :: soilstore_surf ! Initial water stored in soil beneath `Bldgs` surface
      REAL(KIND(1D0)), DIMENSION(nsurf) :: state_surf ! Initial wetness condition on SUEWS land covers.
      REAL(KIND(1D0)), DIMENSION(9) :: WUDay_id ! Daily water use for EveTr, DecTr, Grass [mm]
      REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: soilstore_roof ! Soil moisture of roof [mm]
      REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: state_roof ! wetness status of roof [mm]
      REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: soilstore_wall ! Soil moisture of wall [mm]
      REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: state_wall ! wetness status of wall [mm]
   CONTAINS
      PROCEDURE :: allocHydro => allocHydroState_c
      PROCEDURE :: deallocHydro => deallocHydroState_c
   END TYPE HYDRO_STATE

   TYPE, PUBLIC :: HEAT_STATE
      REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: temp_roof ! interface temperature between depth layers in roof [degC]
      REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: temp_wall ! interface temperature between depth layers in wall [degC]
      REAL(KIND(1D0)), DIMENSION(:, :), ALLOCATABLE :: temp_surf ! interface temperature between depth layers [degC]
      REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: tsfc_roof ! roof surface temperature [degC]
      REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: tsfc_wall ! wall surface temperature [degC]
      REAL(KIND(1D0)), DIMENSION(:), ALLOCATABLE :: tsfc_surf ! surface temperature [degC]
   CONTAINS
      PROCEDURE :: allocHeat => allocHeatState_c
      PROCEDURE :: deallocHeat => deallocHeatState_c
   END TYPE HEAT_STATE

CONTAINS
   SUBROUTINE allocate_hydro_state(self, nlayer)
      IMPLICIT NONE
      TYPE(HYDRO_STATE), INTENT(inout) :: self
      INTEGER, INTENT(in) :: nlayer
      !
      IF (ALLOCATED(self%soilstore_roof)) THEN
         DEALLOCATE (self%soilstore_roof)
      END IF

      IF (ALLOCATED(self%state_roof)) THEN
         DEALLOCATE (self%state_roof)
      END IF

      IF (ALLOCATED(self%soilstore_wall)) THEN
         DEALLOCATE (self%soilstore_wall)
      END IF

      IF (ALLOCATED(self%state_wall)) THEN
         DEALLOCATE (self%state_wall)
      END IF

      ALLOCATE (self%soilstore_roof(nlayer))
      ALLOCATE (self%state_roof(nlayer))
      ALLOCATE (self%soilstore_wall(nlayer))
      ALLOCATE (self%state_wall(nlayer))
      !
   END SUBROUTINE allocate_hydro_state

   SUBROUTINE allocHydroState_c(self, nlayer)
      IMPLICIT NONE
      CLASS(HYDRO_STATE), INTENT(inout) :: self
      INTEGER, INTENT(in) :: nlayer
      !
      CALL allocate_hydro_state(self, nlayer)
      !
   END SUBROUTINE allocHydroState_c

   SUBROUTINE dealloc_hydro_state(self)
      IMPLICIT NONE
      TYPE(HYDRO_STATE), INTENT(inout) :: self
      !
      IF (ALLOCATED(self%soilstore_roof)) THEN
         DEALLOCATE (self%soilstore_roof)
      END IF

      IF (ALLOCATED(self%state_roof)) THEN
         DEALLOCATE (self%state_roof)
      END IF

      IF (ALLOCATED(self%soilstore_wall)) THEN
         DEALLOCATE (self%soilstore_wall)
      END IF

      IF (ALLOCATED(self%state_wall)) THEN
         DEALLOCATE (self%state_wall)
      END IF
      !
   END SUBROUTINE dealloc_hydro_state

   SUBROUTINE deallocHydroState_c(self)
      IMPLICIT NONE
      CLASS(HYDRO_STATE), INTENT(inout) :: self
      !
      CALL dealloc_hydro_state(self)
      !
   END SUBROUTINE deallocHydroState_c

   SUBROUTINE allocate_heat_state(self, num_surf, num_layer, num_depth)
      IMPLICIT NONE
      TYPE(HEAT_STATE), INTENT(inout) :: self
      INTEGER, INTENT(in) :: num_surf, num_layer, num_depth
      !
      IF (ALLOCATED(self%temp_roof)) THEN
         DEALLOCATE (self%temp_roof)
      END IF

      IF (ALLOCATED(self%temp_wall)) THEN
         DEALLOCATE (self%temp_wall)
      END IF

      IF (ALLOCATED(self%tsfc_roof)) THEN
         DEALLOCATE (self%tsfc_roof)
      END IF

      IF (ALLOCATED(self%tsfc_wall)) THEN
         DEALLOCATE (self%tsfc_wall)
      END IF

      IF (ALLOCATED(self%tsfc_surf)) THEN
         DEALLOCATE (self%tsfc_surf)
      END IF

      IF (ALLOCATED(self%temp_surf)) THEN
         DEALLOCATE (self%temp_surf)
      END IF

      ALLOCATE (self%temp_roof(num_layer, num_depth))
      ALLOCATE (self%temp_wall(num_layer, num_depth))
      ALLOCATE (self%tsfc_roof(num_layer))
      ALLOCATE (self%tsfc_wall(num_layer))
      ALLOCATE (self%tsfc_surf(num_surf))
      ALLOCATE (self%temp_surf(num_surf, num_depth))
      !
   END SUBROUTINE allocate_heat_state

   SUBROUTINE allocHeatState_c(self, num_surf, num_layer, num_depth)
      IMPLICIT NONE
      CLASS(HEAT_STATE), INTENT(inout) :: self
      INTEGER, INTENT(in) :: num_surf, num_layer, num_depth
      !
      CALL allocate_heat_state(self, num_surf, num_layer, num_depth)
      !
   END SUBROUTINE allocHeatState_c

   SUBROUTINE deallocHeatState_c(self)
      IMPLICIT NONE
      CLASS(HEAT_STATE), INTENT(inout) :: self
      !
      CALL dealloc_heat_state(self)
      !
   END SUBROUTINE deallocHeatState_c

   SUBROUTINE dealloc_heat_state(self)
      IMPLICIT NONE
      TYPE(HEAT_STATE), INTENT(inout) :: self
      !
      IF (ALLOCATED(self%temp_roof)) THEN
         DEALLOCATE (self%temp_roof)
      END IF

      IF (ALLOCATED(self%temp_wall)) THEN
         DEALLOCATE (self%temp_wall)
      END IF

      IF (ALLOCATED(self%tsfc_roof)) THEN
         DEALLOCATE (self%tsfc_roof)
      END IF

      IF (ALLOCATED(self%tsfc_wall)) THEN
         DEALLOCATE (self%tsfc_wall)
      END IF

      IF (ALLOCATED(self%tsfc_surf)) THEN
         DEALLOCATE (self%tsfc_surf)
      END IF

      IF (ALLOCATED(self%temp_surf)) THEN
         DEALLOCATE (self%temp_surf)
      END IF
      !
   END SUBROUTINE dealloc_heat_state

END MODULE SUEWS_HydroHeat_DTS
