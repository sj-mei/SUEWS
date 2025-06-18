# SUEWS × Claude Code Integration Guide

A comprehensive guide for setting up and using Claude Code with SUEWS for enhanced academic research and development productivity.

## Overview

This integration provides a complete, containerised development environment for SUEWS that includes:

- **Isolated Environment**: Complete SUEWS development stack without affecting your macOS
- **Intelligent Assistance**: Claude Code for code development, debugging, and research
- **Reproducible Research**: Version-controlled environment ensuring consistent results
- **Academic Workflow**: Optimised for research, documentation, and publication
- **Safety First**: Secure credential management and resource isolation

## Quick Start

### 1. One-Time Setup

```bash
# From your SUEWS repository root
./claude-dev/setup-claude-dev.sh
```

This script will:
- Install Claude Code Sandbox
- Create necessary directories
- Set up configuration files
- Create convenience scripts

### 2. Start Development

```bash
./run-claude-dev.sh
```

This launches a complete SUEWS development environment with Claude Code integration.

## Architecture

### Container Environment

The Docker container includes:

- **Ubuntu 22.04** base system
- **GFortran 11** optimised for SUEWS compilation
- **Python 3.12** with complete scientific stack
- **Micromamba** for environment management
- **Development tools**: pytest, ruff, documentation tools
- **SUEWS dependencies**: All required libraries and tools

### File Structure

```
SUEWS/
├── claude-dev/                      # Claude integration files
│   ├── Dockerfile.claude-dev        # Container definition
│   ├── claude-sandbox.config.json   # Sandbox configuration
│   └── setup-claude-dev.sh          # Setup script
├── .env                            # Environment variables (create this)
├── run-claude-dev.sh               # Start development
├── stop-claude-dev.sh              # Stop environment
├── cleanup-claude-dev.sh           # Clean up resources
├── CLAUDE_DEV_GUIDE.md            # Detailed usage guide
├── notebooks/                      # Jupyter notebooks (optional mount)
├── configs/                        # Custom configurations (optional mount)
├── outputs/                        # Simulation outputs (optional mount)
└── external-data/                  # External datasets (optional mount)
```

## Development Workflow

### Inside the Container

Once started, you have access to:

1. **Complete SUEWS environment**: Source code, test data, documentation
2. **Development tools**: make, pytest, jupyter, documentation server
3. **Claude Code assistance**: Intelligent code development and debugging
4. **Git integration**: Version control with automatic branch management

### Common Development Tasks

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
- **Debugging**: Identifying and fixing issues in Fortran and Python code
- **Testing**: Creating comprehensive test suites and interpreting results
- **Documentation**: Writing clear, academic-standard documentation
- **Configuration**: Setting up complex SUEWS configurations
- **Data Analysis**: Analysing simulation results and creating visualisations
- **Research Workflow**: Academic publication and research assistance

## Advanced Features

### Performance Profiling

Built-in tools for performance analysis:

```bash
# Memory profiling
python -m memory_profiler simulation_script.py

# Line-by-line profiling
kernprof -l -v simulation_script.py

# CPU profiling
python -m cProfile -o profile.stats simulation_script.py
```

### Jupyter Integration

Access Jupyter notebooks for analysis:

```bash
# Start Jupyter server
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser

# Access at: http://localhost:8888
```

### Custom Data Integration

Mount external datasets:

```bash
# Place data in external-data/ directory
# Access inside container at /workspace/external-data/
```

### Multi-Site Analysis

Use the extensive test data for multi-site studies:

- **Single-grid data**: `test/data_test/single-grid/`
- **Multi-grid data**: `test/data_test/multi-grid/`
- **Benchmark data**: `test/benchmark1/`

## Configuration Reference

### Sandbox Configuration

The `claude-dev/claude-sandbox.config.json` file controls:

- **Container resources**: CPU, memory, storage limits
- **File mounts**: Host directories accessible in container
- **Environment variables**: Development settings
- **Setup commands**: Automatic environment preparation
- **Security settings**: Isolation and safety protocols

### Environment Variables

Optional variables in `.env`:

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

Default resource limits:

- **CPU**: 6 cores
- **Memory**: 12GB
- **Shared memory**: 2GB
- **Process limit**: 1024

Adjust in `claude-sandbox.config.json` based on your system.

## Academic Research Integration

### Publication Workflow

1. **Development**: Use Claude Code for rapid prototyping and testing
2. **Analysis**: Jupyter notebooks for data analysis and visualisation
3. **Documentation**: Automatic documentation generation with Sphinx
4. **Reproducibility**: Version-controlled environment and configurations
5. **Collaboration**: Git integration for team development

### British Academic Standards

All outputs follow British academic conventions:

- **Language**: British English spelling and grammar
- **Citations**: Academic citation standards
- **Figures**: Publication-ready visualisations
- **Documentation**: Clear, academic-style documentation

### Data Management

- **Test data**: Extensive meteorological datasets included
- **Sample configurations**: Ready-to-use SUEWS setups
- **Output organisation**: Structured output directories
- **Version control**: All changes tracked and versioned

## Security and Safety

### Configuration Management

- **Development settings**: Optional configuration in `.env` file
- **Build configuration**: Parallel job settings and optimisation flags
- **Logging**: Configurable log levels for debugging
- **Access control**: Read-only mounts for sensitive data

### Resource Limits

- **Memory limits**: Prevent runaway processes
- **CPU limits**: Maintain system responsiveness  
- **Process limits**: Control resource usage
- **Network isolation**: Optional network restriction

### Data Protection

- **Read-only mounts**: Protect original datasets
- **Backup integration**: Git version control
- **Cleanup tools**: Safe environment disposal
- **Access logging**: Development activity tracking

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

# Adjust resource limits in claude-sandbox.config.json
```

### Getting Help

1. **Documentation**: `make livehtml` for complete SUEWS documentation
2. **Claude Code**: Ask questions directly to Claude for immediate help
3. **GitHub Issues**: Report bugs and request features
4. **Academic Support**: Contact SUEWS development team

## Best Practices

### Development

1. **Use version control**: Commit changes frequently
2. **Test thoroughly**: Run tests before major changes
3. **Document changes**: Update relevant documentation
4. **Profile performance**: Optimise code with built-in tools
5. **Follow conventions**: Use British English and academic standards

### Research

1. **Reproducible workflows**: Use version-controlled configurations
2. **Data organisation**: Structure outputs for easy analysis
3. **Documentation**: Maintain clear research records
4. **Collaboration**: Use Git for team development
5. **Publication**: Generate publication-ready outputs

### Safety

1. **Secure credentials**: Never commit API keys
2. **Resource monitoring**: Watch memory and CPU usage
3. **Regular cleanup**: Use cleanup scripts to free resources
4. **Backup data**: Version control all important work
5. **Environment isolation**: Keep development isolated from host

## Support and Community

- **SUEWS Repository**: https://github.com/UMEP-dev/SUEWS
- **Documentation**: https://suews.readthedocs.io/
- **Issues and Support**: GitHub Issues
- **Academic Community**: SUEWS research network
- **Claude Code**: https://claude.ai/code

---

*This integration guide is part of the SUEWS project and follows British academic standards. For technical support, please use the GitHub issue tracker or contact the development team.*