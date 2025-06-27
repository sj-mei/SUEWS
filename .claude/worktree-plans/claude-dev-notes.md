# SUEWS Development Plans Overview

This document provides a summary of all active development branches and their associated GitHub issues.

## Active Branches and Issues

### 1. feature/core-runtime-fixes
**Focus**: Critical runtime bugs and stability issues  
**Key Issues**:
- #391: Test reinitialisation of SuPy (STEBBS)
- #406: LAI and SMD not correct using sample data
- #335: Water balance issue (verify fix)
- #354: STEBBS water system issues (verify fix)

### 2. feature/adjust-default-values
**Focus**: Removing problematic default values  
**Key Issues**:
- #428: Removing or adjusting default values (PRIMARY)
- #412: Update sample data with benchmark

### 3. feature/infrastructure-enhancements
**Focus**: Build system, testing, and infrastructure improvements  
**Key Issues**:
- #400: Pydantic Hierarchy
- #392: Pydantic-related cleanup
- #385: New Pytest benchmarking
- #350: Couple SPARTACUS to STEBBS
- #360: Enable multiple building archetypes

### 4. feature/fast-dev-build
**Focus**: Fast builds and documentation overhaul  
**Key Issues**:
- #411: Documentation system revision (active work)
- #343: Update docs site
- #365: Interactive config editor
- #328: v2025 release checklist

### 5. feature/rsl-physics-fixes
**Focus**: Roughness Sublayer physics corrections  
**Key Issues**:
- #419: MOST in RSL for short buildings (PRIMARY)
- #338: High z0 for low FAI cases
- #349: RSL-STEBBS coupling (verify)

### 6. feature/supy-data-processing
**Focus**: Data processing and I/O improvements  
**Key Issues**:
- #408: Full output from save_supy
- #417: Output building OHM coefficients
- #412: Update sample data
- #406: LAI/SMD preprocessing
- #353: Unit standardisation

### 7. feature/hide-internal-options
**Focus**: Hide developer options from user docs  
**Key Issues**: 
- No specific issue number (initiated from discussion)
- Already has detailed plan in place

## Priority Order

1. **High Priority**: 
   - core-runtime-fixes (stability issues)
   - adjust-default-values (breaking changes)
   - fast-dev-build (v2025 release)

2. **Medium Priority**:
   - rsl-physics-fixes (physics accuracy)
   - supy-data-processing (user experience)

3. **Lower Priority**:
   - infrastructure-enhancements (long-term improvements)
   - hide-internal-options (documentation cleanup)

## Coordination Notes

- Several branches touch the pydantic data model - coordinate changes
- Documentation updates needed across multiple branches
- Test coverage should be added with each fix
- Consider release timing for breaking changes