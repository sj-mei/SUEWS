# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SUEWS (Surface Urban Energy and Water balance Scheme) is a physics-based model for simulating urban surface energy and water balance. It combines a Fortran physics engine with a Python interface called `supy`.

## Build System and Commands

### Essential Development Commands

```bash
# Setup development environment
mamba env create -f env.yml
mamba activate suews-dev

# Initialize SPARTACUS submodule (required before first build)
git submodule init
git submodule update

# Development build (fast, no tests)
make dev

# Full build with tests
make

# Run tests only
make test

# Build documentation
make docs

# Live preview documentation with auto-rebuild
make livehtml

# Generate configuration schema from Pydantic models
make schema

# Process CSV files for documentation
make proc-csv

# Start SUEWS Configuration UI server (http://localhost:8080)
make config-ui

# Convert legacy table-based input to YAML
suews-convert to-yaml -i path/to/input/dir -o config.yml

# Convert with version upgrade (e.g., from 2020a)
suews-convert to-yaml -i path/to/input/dir -o config.yml -f 2020a

# Clean all build artifacts  
make clean

# Build distribution wheels
make wheel
```

### Key Build System Details

- **Build system**: Meson with meson-python backend
- **Languages**: Fortran 90 + C + Python
- **Primary compiler**: gfortran ≥9.3.0 (Intel ifort also supported)
- **Development install**: Uses `--no-build-isolation --editable` flags for faster iteration

## Architecture

### Three-Component Structure

1. **`src/suews/`**: Fortran physics engine
   - Modular design with separate modules for different physics:
     - `suews_phys_*`: Core physics (evaporation, snow, radiation, etc.)
     - `suews_ctrl_*`: Control flow and I/O handling  
     - `suews_util_*`: Utilities (datetime, meteorology, etc.)
   - External dependency: SPARTACUS radiative transfer model (git submodule)
   - Builds static libraries: `libsuewsdriver.a`, `libsuewsphys.a`, `libsuewsutil.a`

2. **`src/supy/`**: Python interface and utilities
   - Modern Pydantic-based configuration in `data_model/` subpackage
   - Command-line tools: `suews-run`, `suews-convert`, `to_yaml`
   - Utility modules in `util/`: ERA5 processing, plotting, gap filling, TMY data, etc.
   - Sample configurations in `sample_run/`

3. **`src/supy_driver/`**: F2Py wrapper connecting Fortran to Python
   - Uses f90wrap for automated wrapper generation
   - Contains build scripts and patches for the F2Py interface

### Data Model Transition

The project recently migrated from `_config.py` to a new `data_model/` subpackage with Pydantic-based classes. When working with configuration classes, use the new `data_model` module.

## Testing

### Test Framework
- **Framework**: pytest with coverage reporting
- **Location**: `/test/` directory
- **Command**: `make test` or `python -m pytest test -v`

### Test Types
- Single-grid multi-year simulations
- Multi-grid multi-year simulations  
- Physics scheme validation
- Regression testing against baseline results (pickled reference data)

### Test Data
- Extensive meteorological datasets (ERA5, TMY) in `test/data_test/`
- Benchmark data for regression tests in `test/benchmark1/`
- Reference results stored as pickled DataFrames

## Development Guidelines

### Branch Strategy
- **Main branch**: `master` (stable, protected)
- **Development**: Feature branches with pull requests
- **Current**: `sunt05/issue343` (active development branch)

### Branch Development Documentation

**CRITICAL**: When working on feature branches, maintain proper documentation throughout the development cycle:

**Branch Naming Conventions:**
- **Feature branches**: `feature/description` (e.g., `feature/suews-yaml-config`)
- **Bug fixes**: `fix/description` (e.g., `fix/docker-permissions`)
- **Documentation**: `docs/description` (e.g., `docs/api-reference`)
- **Refactoring**: `refactor/description` (e.g., `refactor/makefile-cleanup`)
- **Testing**: `test/description` (e.g., `test/regression-suite`)
- **Chore/maintenance**: `chore/description` (e.g., `chore/dependency-update`)

**Documentation Requirements for Branch Work:**

1. **Branch Creation**:
   - Document the purpose and scope in the first commit message
   - Update CHANGELOG.md with `[Unreleased]` entries immediately
   - Create or reference the related GitHub issue

2. **During Development**:
   - Keep CHANGELOG.md updated with each significant change
   - Use descriptive commit messages following conventional format:
     ```
     type(scope): subject
     
     body (optional)
     
     footer (optional)
     ```
   - Examples:
     - `feat(supy): add YAML configuration support`
     - `fix(build): resolve meson conflicts in virtual environments`
     - `docs(api): update docstrings for data_model classes`

3. **Before Pull Request**:
   - Ensure all changes are documented in CHANGELOG.md
   - Update relevant documentation (README, API docs, etc.)
   - Run tests and document any new test cases added
   - Create comprehensive PR description summarising:
     - What changed and why
     - How to test the changes
     - Any breaking changes or migration steps
     - Related issues/PRs

4. **Branch Progress Tracking**:
   - For long-running branches, maintain a branch-specific log:
     ```markdown
     ## Branch: feature/suews-yaml-config
     
     ### Session 1 (2025-06-21)
     - Implemented basic YAML parser
     - Added Pydantic models for configuration
     - Updated build system to handle new dependencies
     
     ### Session 2 (2025-06-22)
     - Added conversion utilities from legacy format
     - Created comprehensive test suite
     - Fixed edge cases in array handling
     ```

5. **Merge Documentation**:
   - Update CHANGELOG.md to consolidate all branch changes
   - Ensure version numbers are updated if applicable
   - Archive branch progress logs for future reference

### Code Quality
- **Fortran**: Use fprettify for automatic formatting
- **Python**: Follow standard conventions
- **Tests**: Always run `make test` before submitting changes
- **Coverage**: pytest-cov provides coverage reporting

### Version Management
- Git-based versioning via `get_ver_git.py`
- Dynamic version in `pyproject.toml`
- Version info embedded in Fortran code

## Key File Locations

### Configuration
- **Python package**: `pyproject.toml` (dependencies, scripts, metadata)
- **Environment**: `env.yml` (unified conda environment for development and documentation)
- **Build**: `meson.build` (primary build configuration)
- **Legacy**: `Makefile` (convenience wrapper)

### Input/Output
- **Sample inputs**: `src/supy/sample_run/`
- **Schema definitions**: `schema/` directory with YAML input specifications
- **File format rules**: `src/supy/util/rules.csv`

### Documentation
- **Source**: `docs/source/` (Sphinx-based with reStructuredText and Markdown)
- **Build**: `make docs` generates HTML, PDF, and API documentation
- **Live preview**: `make livehtml` with auto-rebuild during development
- **Configuration**: `docs/source/conf.py` with extensive customizations
- **See**: `docs/README.md` for complete documentation system overview

## Platform Support

### Supported Platforms
- Linux (x86_64)
- macOS (x64 and ARM64)
- Windows (AMD64)
- Python 3.10-3.13 (3.12 recommended)

### CI/CD
- Multi-platform wheel building via GitHub Actions
- Automated testing on all supported platforms
- PyPI publishing on tagged releases
- Code formatting checks with fprettify

## Common Development Patterns

### Adding New Physics Modules
1. Create `suews_phys_newmodule.f95` in `src/suews/src/`
2. Add corresponding module interface in Fortran
3. Update `src/supy/` for Python interface if needed
4. Add tests in `test/`

### Adding Utility Functions
- **Fortran utilities**: Add to appropriate `suews_util_*.f95`
- **Python utilities**: Add to `src/supy/util/`
- **Data processing**: Consider ERA5, TMY, or atmospheric modules

### Working with Configuration
- Use classes from `src/supy/data_model/` for type-safe configuration
- YAML schema definitions are in `schema/input/`
- Sample configurations demonstrate proper usage patterns
- **Legacy Migration**: Use `suews-convert to-yaml` to convert old table-based inputs
- **Interactive Builder**: Use web UI at `docs/source/_static/index.html` for guided configuration

## Documentation System

### Documentation Architecture
The `docs/` directory contains a sophisticated Sphinx-based documentation system:

- **Build system**: Sphinx with book theme and 20+ extensions
- **Content format**: reStructuredText (.rst) and Markdown (.md)
- **Multi-format output**: HTML (web), PDF (LaTeX), API docs (Doxygen)
- **Live development**: `make livehtml` for real-time preview

### Documentation Structure
```
docs/
├── source/                 # Main content
│   ├── conf.py            # Sphinx configuration (900+ lines)
│   ├── index.rst          # Documentation entry point
│   ├── input_files/       # Input format specifications
│   ├── output_files/      # Output format specifications
│   ├── contributing/      # Development guides
│   ├── version-history/   # Release notes
│   └── images/           # 300+ documentation images
├── build/                 # Generated output
├── suews-config-ui/      # React configuration interface
└── README.md             # Complete system overview
```

### Essential Documentation Commands
```bash
# From docs/ directory
make html              # Build HTML documentation
make livehtml         # Live preview with auto-rebuild
make latexpdf         # Generate PDF manual
make proc-csv         # Process parameter tables
```

### Adding Documentation
1. **New content**: Place RST/MD files in appropriate `source/` subdirectories
2. **Images**: Add to `source/images/` with descriptive names
3. **References**: Update `source/assets/refs/` BibTeX files
4. **Tables**: Use CSV + RST pattern for parameter documentation
5. **API docs**: Doxygen auto-generates from Fortran/C code comments

### RST Syntax Guidelines
**CRITICAL**: All `.rst` files MUST use reStructuredText syntax, NOT Markdown!

- **Headings**: Use `=`, `-`, `~`, `^` underlines (not `#` markers)
- **Links**: Use `:doc:`, `:ref:`, or `link text <url>`_ (not `[text](url)`)
- **Code blocks**: Use `.. code-block:: language` (not triple backticks)
- **Lists**: Use `-` or `*` for bullets, numbers for ordered (same as Markdown)
- **Emphasis**: Use `*italic*` and `**bold**` (same as Markdown)
- **Directives**: Use `.. directive::` format (e.g., `.. note::`, `.. warning::`)
- **Cross-references**: Use `:option:`, `:doc:`, `:ref:` for internal links

### Language and Writing Style
**CRITICAL**: Always use British English for all documentation, annotations, comments, and text content!

**British English Requirements:**
- **Spelling**: Use British spellings throughout:
  - "modelling" not "modeling"  
  - "colour" not "color"
  - "realise" not "realize"
  - "optimise" not "optimize"  
  - "behaviour" not "behavior"
  - "centre" not "center"
  - "analyse" not "analyze"
  - "characterise" not "characterize"
  - "parameterisation" not "parameterization"
  - "visualisation" not "visualization"
- **Grammar**: Use British conventions for punctuation and phrasing
- **Consistency**: Apply British English across all documentation, code comments, docstrings, and user-facing text
- **Technical terms**: Use British scientific terminology where applicable

### Documentation Features
- **Academic citations**: Custom bibliography with author-year format
- **Cross-references**: Auto-linking to pandas, numpy, xarray documentation
- **Interactive elements**: Web configuration UI integration
- **Multi-language**: Support for English with internationalization framework
- **GitHub integration**: Issue templates and edit-on-GitHub functionality
- **Search**: Full-text search with highlighting

### Documentation Dependencies
All specified in the unified `env.yml`:
- **Core**: sphinx, sphinx-book-theme, sphinx-autobuild
- **Extensions**: bibtex, panels, autodoc, napoleon
- **Processing**: doxygen, pandoc, node.js
- **Academic**: pybtex for bibliography management

See `docs/README.md` for comprehensive documentation system details.

## Configuration UI System

### Schema-Driven Configuration Builder

SUEWS includes a web-based configuration builder that generates YAML configurations from the Pydantic data model:

- **Location**: `docs/source/_static/` contains the configuration UI
- **Schema Generation**: `docs/gen_schema.py` extracts JSON Schema from SUEWSConfig
- **Interface Files**:
  - `config-builder.html`: Main configuration interface
  - `config-builder.js`: Interactive form logic with validation  
  - `config-builder.css`: UI styling and responsive design
  - `index.html`: Landing page with feature overview
  - `suews-config-schema.json`: JSON Schema (auto-generated)

### Key Features

- **Real-time Validation**: Live validation against Pydantic schema
- **Export Formats**: Generate YAML or JSON files ready for SUEWS
- **Bootstrap UI**: Modern, responsive interface with tabbed organization
- **Reference Support**: Handles ValueWithDOI fields for scientific references
- **Array Management**: Dynamic addition/removal of sites and profiles

### Development Workflow

1. **Modify Data Model**: Update Pydantic classes in `src/supy/data_model/`
2. **Build SUEWS**: Run `make clean && make dev` to rebuild supy
3. **Generate Schema**: Run `cd docs && python gen_schema.py`
4. **Test UI**: Open `docs/source/_static/index.html` to test configuration builder

### Schema Generation Details

- **Single Source**: Only `docs/source/_static/suews-config-schema.json` is created
- **Automatic**: Schema reflects current Pydantic model definitions
- **Validation**: Full JSON Schema with types, constraints, and descriptions
- **Integration**: Used by both config builder and documentation

## Changelog Management

### Automatic Changelog Updates

**CRITICAL**: When making significant changes to SUEWS, automatically update `CHANGELOG.md` following these guidelines:

**What constitutes significant changes:**
- New features or functionality
- Breaking changes or API modifications
- Bug fixes that affect user workflows
- Performance improvements
- Documentation restructuring
- Build system modifications
- New command-line tools or options
- Configuration schema changes

**Changelog format to follow:**
```markdown
## [Unreleased]

### Added
- New features and functionality

### Changed
- Changes to existing functionality

### Deprecated
- Features marked for removal in future versions

### Removed
- Features removed in this version

### Fixed
- Bug fixes

### Security
- Security-related changes
```

**When to update:**
- **During development**: Add entries to `[Unreleased]` section immediately when implementing changes
- **After completing multi-step changes**: Summarise the complete feature/improvement
- **Before pull requests**: Ensure all changes are documented
- **Never retroactively**: Always update changelog as part of the change implementation

**Writing style:**
- Use British English spellings (realise, optimise, colour, etc.)
- Write from user perspective ("Added support for..." not "We added...")
- Be specific and actionable
- Include relevant command examples where helpful
- Reference issue numbers when applicable

**Example entry for recent work:**
```markdown
## [Unreleased]

### Added
- Cross-platform isolated build directories to prevent environment conflicts
- Automatic environment detection in `make dev` command
- New Makefile target: `make deactivate`
- Enhanced help system with `make help` showing Quick Start guide

### Changed
- Improved Makefile with unified development workflow
- Enhanced cross-platform compatibility for Windows, macOS, and Linux
- Optimised build process with isolated build directories in `/tmp/suews-builds/`

### Fixed
- Resolved meson build conflicts between different Python environments
- Fixed numpy path issues when using virtual environments within project directory
```

**Integration with development workflow:**
1. **Before implementing changes**: Consider if changelog entry will be needed
2. **During implementation**: Add/update changelog entries as you work
3. **Before committing**: Verify changelog reflects all changes made
4. **During PR review**: Ensure changelog entries are accurate and complete

## Progress Logging Requirements

### Automatic Session Documentation

**CRITICAL**: At the end of significant development sessions, automatically create a detailed progress log and update relevant changelog files:

**What to document:**
- **Major features implemented** (new commands, build system changes, environment setup)
- **Bug fixes and technical solutions** (build conflicts, path issues, cross-platform fixes)
- **User experience improvements** (help systems, convenience commands, workflow enhancements)
- **Repository maintenance** (gitignore updates, documentation improvements)
- **Technical architecture changes** (isolated builds, environment detection, worktree support)

**Where to document:**
1. **`CHANGELOG.md`**: User-facing changes following the established format
2. **Progress log in session**: Comprehensive technical summary for development records
3. **`CLAUDE.md` updates**: Any new development patterns or guidelines discovered

**Progress log format:**
```markdown
## Progress Log - [Session Topic]

### ✅ **[Category Name]**
- Specific accomplishment with technical details
- Impact and benefits achieved
- Cross-references to related changes

### 🎯 **Key Achievements**
- High-level summary of major accomplishments
- User benefits and workflow improvements

### 📊 **New Features/Commands Added**
- Specific commands, targets, or functionality

### 🔧 **Technical Solutions Implemented**
- Detailed technical fixes and their rationale
- Architecture decisions and implementation details
```

**Trigger criteria for progress logging:**
- Multi-step feature implementations
- Build system modifications
- New development workflows or tools
- Cross-platform compatibility improvements
- Major bug fixes affecting user workflows
- Repository structure or documentation overhauls

**Files to check for updates:**
- `CHANGELOG.md` (primary user-facing changes)
- `CLAUDE.md` (development guidelines and patterns)
- `README.md` (if user-facing features change)
- `docs/` files (if documentation workflow changes)
- Any project-specific changelog or history files

**Timing:**
- **End of major development sessions**: Create comprehensive progress log
- **Before completing significant features**: Ensure all changes are documented
- **When switching contexts**: Document progress before moving to different work