# SUEWS Claude Code Development Environment

A comprehensive containerised development environment for SUEWS with Claude Code integration, providing intelligent AI-assisted development capabilities for enhanced productivity.

## Overview

This integration provides:

- **Isolated Environment**: Complete SUEWS development stack without affecting your macOS
- **Intelligent Assistance**: Claude Code for code development, debugging, and research
- **Reproducible Research**: Version-controlled environment ensuring consistent results
- **Academic Workflow**: Optimised for research, documentation, and publication
- **British Standards**: All outputs follow British academic conventions

## Quick Start

### From Repository Root

```bash
# Using Makefile (recommended) - creates workspace outside Dropbox
make claude-dev

# Using custom location (e.g., for specific project folder)
make claude-dev LOCATION=/tmp/suews-workspace

# Or using setup script directly (from current location)
./claude-dev/setup-claude-dev.sh
```

### Start Development

After running `make claude-dev`, navigate to the created workspace:

```bash
# Default location
cd ~/claude-suews-workspace/SUEWS && ./run-claude-dev.sh

# Or use the path shown in setup output
cd /path/to/your/workspace/SUEWS && ./run-claude-dev.sh
```

### Dropbox Compatibility

**Important**: Claude Code Sandbox cannot run properly from within Dropbox folders due to file monitoring conflicts. The `make claude-dev` recipe automatically:

1. Creates a workspace outside Dropbox (default: `~/claude-suews-workspace`)
2. Clones the SUEWS repository to that location
3. Sets up Claude Code environment there

**Recommended workflow for Dropbox users**:
- Use `make claude-dev` instead of running setup directly
- Work in the cloned workspace, not the Dropbox folder
- Push changes back to the main repository when ready

### Workspace Management

The `make claude-dev` recipe handles workspace management intelligently:

- **First run**: Clones SUEWS repository and sets up environment
- **Subsequent runs**: Updates existing repository and refreshes setup
- **Branch handling**: Automatically switches to `claude-dev` branch if available
- **Custom locations**: Use `LOCATION=/your/path` to specify workspace directory

**Example workflows**:
```bash
# Set up in default location
make claude-dev

# Set up in temporary directory
make claude-dev LOCATION=/tmp/suews-dev

# Set up in project-specific folder
make claude-dev LOCATION=/Users/username/Projects/suews-research
```

## Files in This Directory

- **`Dockerfile.claude-dev`**: Production-ready container with complete SUEWS environment
- **`claude-sandbox.config.json`**: Configuration for Claude Code Sandbox
- **`setup-claude-dev.sh`**: Automated setup script
- **`.dockerignore`**: Optimised Docker build context
- **`README.md`**: This comprehensive guide

## Container Environment

### Base System
- **OS**: Ubuntu 22.04 LTS
- **Compiler**: GFortran 11 (optimised for SUEWS)
- **Python**: 3.12 with complete scientific stack
- **Package Manager**: Micromamba for fast, reliable environments

### Development Tools
- **Build System**: Meson with ninja backend
- **Testing**: pytest with coverage reporting
- **Code Quality**: ruff, black, mypy
- **Documentation**: Complete Sphinx ecosystem
- **Profiling**: memory_profiler, line_profiler, cProfile

### SUEWS Integration
- **Source Code**: Complete repository with submodules
- **Test Data**: All benchmark and sample data included
- **Sample Configurations**: Ready-to-use SUEWS setups
- **Development Build**: Optimised for rapid iteration

### Claude Code Features
- **Intelligent Assistance**: AI-powered code development and debugging
- **Git Integration**: Automatic branch management
- **Documentation**: Automated generation and maintenance
- **Testing**: Intelligent test creation and interpretation
- **Research Support**: Academic publication workflow

## Development Workflow

### Inside the Container

Once started, you have access to:

1. **Complete SUEWS environment**: Source code, test data, documentation
2. **Development tools**: make, pytest, jupyter, documentation server
3. **Claude Code assistance**: Intelligent development support
4. **Git integration**: Version control with automatic branch management

### Common Commands

#### Building SUEWS

```bash
# Quick development build
make dev

# Full build with tests
make

# Clean build
make clean
```

#### Testing

```bash
# Run all tests
make test

# Run specific test
python -m pytest test/test_data_model.py -v

# Run with coverage
python -m pytest test/ --cov=supy
```

#### Documentation

```bash
# Build documentation with live reload
make livehtml

# Build static documentation
make docs

# Process CSV files for documentation
make proc-csv
```

#### Running Simulations

```bash
# Run sample simulation
cd src/supy/sample_run
suews-run .

# Convert legacy configuration
suews-convert to-yaml -i old_config_dir -o new_config.yml

# Run with custom configuration
suews-run path/to/config.yml
```

### Claude Code Assistance

Claude Code provides intelligent assistance for:

- **Code Development**: Writing, refactoring, and optimising SUEWS code
- **Debugging**: Identifying and fixing issues in Fortran and Python
- **Testing**: Creating comprehensive test suites
- **Documentation**: Writing clear, academic-standard documentation
- **Configuration**: Setting up complex SUEWS configurations
- **Data Analysis**: Analysing results and creating visualisations
- **Research Workflow**: Academic publication support

## Configuration

### Sandbox Configuration

The `claude-sandbox.config.json` file controls:

- **Container resources**: CPU, memory, storage limits
- **File mounts**: Host directories accessible in container
- **Environment variables**: Development settings
- **Setup commands**: Automatic environment preparation
- **Security settings**: Isolation and safety protocols

### Environment Variables

Optional configuration in `.env`:

```bash
# Development Configuration
SUEWS_DEV_MODE=true                   # Development features
SUEWS_BUILD_JOBS=4                    # Parallel build jobs
SUEWS_LOG_LEVEL=INFO                  # Logging level

# Optional External Services
# CDSAPI_URL=https://cds.climate.copernicus.eu/api/v2
# CDSAPI_KEY=your-cds-api-key
```

### Resource Configuration

Default limits (adjustable in `claude-sandbox.config.json`):

- **CPU**: 6 cores
- **Memory**: 12GB
- **Shared memory**: 2GB
- **Process limit**: 1024

## Advanced Features

### Performance Profiling

```bash
# Memory profiling
python -m memory_profiler simulation_script.py

# Line-by-line profiling
kernprof -l -v simulation_script.py

# CPU profiling
python -m cProfile -o profile.stats simulation_script.py
```

### Jupyter Integration

```bash
# Start Jupyter server
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser

# Access at: http://localhost:8888
```

### Custom Data Integration

Mount external datasets by placing them in:
- `notebooks/` - Jupyter notebooks
- `configs/` - Custom SUEWS configurations
- `outputs/` - Simulation results
- `external-data/` - External datasets (read-only)

## Academic Research Integration

### Publication Workflow

1. **Development**: Rapid prototyping with Claude Code
2. **Analysis**: Jupyter notebooks for data analysis
3. **Documentation**: Automatic documentation generation
4. **Reproducibility**: Version-controlled environment
5. **Collaboration**: Git integration for teams

### British Academic Standards

- **Language**: British English spelling and grammar
- **Citations**: Academic citation standards
- **Figures**: Publication-ready visualisations
- **Documentation**: Clear, academic-style writing

### Data Management

- **Test data**: Extensive meteorological datasets included
- **Sample configurations**: Ready-to-use setups
- **Output organisation**: Structured directories
- **Version control**: All changes tracked

## Security and Safety

### Configuration Management

- **Development settings**: Optional `.env` configuration
- **Build configuration**: Parallel job settings
- **Logging**: Configurable log levels
- **Access control**: Read-only mounts for data

### Resource Protection

- **Memory limits**: Prevent system exhaustion
- **Process limits**: Control concurrent operations
- **Network isolation**: Optional restriction
- **File system**: Isolated container environment

## Troubleshooting

### Common Issues

#### Container Won't Start

```bash
# Check Docker daemon
docker version

# Check Claude Code Sandbox
claude-sandbox --version

# Verify configuration
cat claude-dev/claude-sandbox.config.json
```

#### Build Failures

```bash
# Check compiler versions
gcc --version
gfortran --version

# Verify dependencies
micromamba list -n suews-dev

# Clean and rebuild
make clean && make dev
```

#### Import Errors

```bash
# Check Python environment
python -c "import supy; print(supy.__version__)"

# Verify installation
pip show supy

# Reinstall in development mode
pip install -e .
```

#### Memory Issues

```bash
# Check available memory
free -h

# Monitor usage
htop

# Adjust limits in claude-sandbox.config.json
```

### Getting Help

1. **Documentation**: `make livehtml` for complete SUEWS docs
2. **Claude Code**: Ask questions directly for immediate help
3. **GitHub Issues**: Report bugs and request features
4. **Academic Support**: Contact SUEWS development team

## Best Practices

### Development

1. **Version control**: Commit changes frequently
2. **Test thoroughly**: Run tests before major changes
3. **Document changes**: Update relevant documentation
4. **Profile performance**: Optimise with built-in tools
5. **Follow conventions**: British English and academic standards

### Research

1. **Reproducible workflows**: Version-controlled configurations
2. **Data organisation**: Structure outputs for analysis
3. **Documentation**: Maintain clear research records
4. **Collaboration**: Use Git for team development
5. **Publication**: Generate publication-ready outputs

### Maintenance

1. **Regular cleanup**: Use cleanup scripts
2. **Monitor resources**: Watch memory and CPU usage
3. **Update dependencies**: Keep environment current
4. **Backup work**: Version control all changes
5. **Environment isolation**: Keep development separate

## Maintenance

### Regular Cleanup

```bash
# Stop all containers
./stop-claude-dev.sh

# Clean up resources
./cleanup-claude-dev.sh

# Remove unused Docker images
docker image prune -f
```

### Updates

```bash
# Update Claude Code Sandbox
npm update -g @anthropic-ai/claude-code @textcortex/claude-code-sandbox

# Rebuild container with latest dependencies
docker build -t suews-claude-dev -f claude-dev/Dockerfile.claude-dev .
```

## Support

- **SUEWS Documentation**: https://suews.readthedocs.io/
- **Claude Code**: https://claude.ai/code
- **GitHub Issues**: https://github.com/UMEP-dev/SUEWS/issues
- **Academic Community**: SUEWS research network

---

*This integration is optimised for academic research and follows British academic standards. For technical support, please use the GitHub issue tracker.*