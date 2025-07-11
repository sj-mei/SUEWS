# Feature: Cross-Platform Deterministic Build

## Context
Implement deterministic build mode for SUEWS to guarantee identical results across different platforms, compilers, and architectures with the same input. This addresses the need for reproducible scientific computations and consistent validation across development environments.

## GitHub Issues
- None directly related (new feature)

## Progress Tracking

### Phase 1: Build Infrastructure
- [x] Create deterministic math module with Kahan summation
- [x] Add deterministic build flags to meson.build
- [x] Add deterministic build flags to Makefile.gfortran
- [x] Set up dedicated worktree with uv environment
- [ ] Test basic compilation with deterministic flags

### Phase 2: Code Modifications
- [ ] Replace SUM() with det_sum() in critical modules
- [ ] Standardize numerical tolerances
- [ ] Replace transcendental functions with portable versions
- [ ] Fix random number seeding
- [ ] Handle big-endian conversion consistently

### Phase 3: Testing Infrastructure
- [ ] Create cross-platform test cases
- [ ] Implement numerical fingerprinting
- [ ] Add CI matrix for Linux/macOS/Windows
- [ ] Create Docker test environment
- [ ] Document performance impact

### Phase 4: Integration
- [ ] Update cibuildwheel configuration
- [ ] Add deterministic mode to documentation
- [ ] Create user guide for reproducible runs
- [ ] Add performance benchmarks

## Key Decisions
- Use Kahan summation for all critical accumulations
- Standardize on -O2 optimization for deterministic mode
- Use statistical fingerprinting instead of bit-for-bit comparison
- Make deterministic mode optional to preserve performance

## Implementation Notes
- Deterministic module created at `src/suews/suews_util_deterministic.f95`
- Compiler flags disable FP contraction and fast math
- Need to handle big-endian conversion flag impact
- Windows ARM testing deferred to later phase

## Files to Modify
- `src/suews/suews_util_deterministic.f95` (created)
- `meson_options.txt` (modified)
- `meson.build` (modified) 
- `src/supy_driver/meson.build` (modified)
- `src/suews/Makefile` (modified)
- `src/suews/Makefile.gfortran` (modified)
- Physics modules for SUM replacement (pending)
- Test suite additions (pending)

## Current Status
Working in dedicated worktree `worktrees/deterministic-build` with uv environment set up. Ready to test compilation and begin code modifications.