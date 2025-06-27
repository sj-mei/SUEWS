# Feature: Infrastructure Enhancements

## Context
This branch focuses on improving the overall infrastructure of SUEWS/SuPy including build system, testing framework, performance optimisations, and development workflow improvements.

## GitHub Issues to Address
- **#400**: Pydantic Hierarchy (enhancement, pydantic)
- **#392**: Pydantic-related cleanup (dev, pydantic)
- **#385**: New Pytest benchmarking (benchmark)
- **#365**: Add an interactive editor for SUEWS config - CLOSED but verify
- **#350**: Couple SPARTACUS to STEBBS (SPARTACUS, STEBBS)
- **#360**: Enable multiple building archetypes (STEBBS)

## Progress Tracking
- [ ] Pydantic infrastructure improvements
  - [ ] Implement proper hierarchy for pydantic models (#400)
  - [ ] Clean up pydantic-related code (#392)
  - [ ] Improve validation performance
  - [ ] Enhance error messages
- [ ] Testing infrastructure
  - [ ] Implement new pytest benchmarking system (#385)
  - [ ] Add performance regression tests
  - [ ] Improve test coverage reporting
  - [ ] Add continuous benchmarking
- [ ] Build system enhancements
  - [ ] Optimise Fortran compilation
  - [ ] Improve Python packaging
  - [ ] Add development mode builds
  - [ ] Enhance CI/CD pipeline
- [ ] Model coupling infrastructure
  - [ ] Prepare infrastructure for SPARTACUS-STEBBS coupling (#350)
  - [ ] Enable multiple building archetypes (#360)
  - [ ] Improve inter-model communication

## Key Decisions
- Maintain backward compatibility for public APIs
- Focus on developer experience improvements
- Prioritise performance for common operations
- Keep infrastructure changes modular and testable

## Implementation Notes
- Pydantic v2 migration should be considered
- Benchmark suite should track key metrics over time
- Build system should support both development and production modes
- Consider using modern Python packaging standards (PEP 517/518)

## Files to Modify
- `src/supy/data_model/` - Pydantic hierarchy improvements
- `test/` - New benchmarking infrastructure
- `pyproject.toml` - Modern packaging configuration
- `meson.build` - Build system enhancements
- `.github/workflows/` - CI/CD improvements
- `src/suews/src/` - Model coupling preparations

## Infrastructure Goals
1. Reduce build times for development
2. Improve error messages and debugging
3. Enable easier testing and benchmarking
4. Support advanced model configurations
5. Streamline release process

## Performance Targets
- 20% reduction in build time
- 10% improvement in validation performance
- Comprehensive benchmark coverage
- Sub-second test suite for unit tests