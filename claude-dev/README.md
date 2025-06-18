# SUEWS Claude Code Docker Integration

This directory contains Docker configuration files for running SUEWS with Claude Code in a containerised development environment.

## Files

- **`Dockerfile.claude-dev`**: Production-ready Docker container with complete SUEWS development environment
- **`claude-sandbox.config.json`**: Configuration for Claude Code Sandbox integration
- **`setup-claude-dev.sh`**: Automated setup script for the entire development environment
- **`.dockerignore`**: Optimised Docker build context exclusions
- **`README.md`**: This documentation file

## Quick Start

1. **Run setup** (from SUEWS repository root):
   ```bash
   ./claude-dev/setup-claude-dev.sh
   ```

2. **Start development**:
   ```bash
   ./run-claude-dev.sh
   ```

## Container Features

### Base Environment
- **OS**: Ubuntu 22.04 LTS
- **Compiler**: GFortran 11 (optimised for SUEWS)
- **Python**: 3.12 with scientific computing stack
- **Package Manager**: Micromamba for fast, reliable package management

### Development Tools
- **Build System**: Meson with ninja backend
- **Testing**: pytest with coverage reporting
- **Code Quality**: ruff, black, mypy
- **Documentation**: Complete Sphinx ecosystem
- **Profiling**: memory_profiler, line_profiler, cProfile

### SUEWS Integration
- **Source Code**: Complete SUEWS repository with submodules
- **Test Data**: All benchmark and sample data included
- **Sample Configurations**: Ready-to-use SUEWS setups
- **Development Build**: Optimised for rapid iteration

### Claude Code Features
- **Intelligent Assistance**: AI-powered code development and debugging
- **Git Integration**: Automatic branch management and commit assistance
- **Documentation**: Automated documentation generation and maintenance
- **Testing**: Intelligent test creation and result interpretation

## Configuration

### Resource Limits
Default container limits (adjustable in `claude-sandbox.config.json`):
- **CPU**: 6 cores
- **Memory**: 12GB
- **Shared Memory**: 2GB
- **Process Limit**: 1024

### File Mounts
Optional directory mounts for external resources:
- `notebooks/` → `/workspace/notebooks` (Jupyter notebooks)
- `configs/` → `/workspace/configs` (Custom SUEWS configurations)
- `outputs/` → `/workspace/outputs` (Simulation results)
- `external-data/` → `/workspace/external-data` (External datasets, read-only)

### Environment Variables
Optional configuration options:
- `SUEWS_DEV_MODE`: Enable development features
- `SUEWS_BUILD_JOBS`: Number of parallel build jobs
- `SUEWS_LOG_LEVEL`: Logging level for debugging

## Development Workflow

### Inside Container

```bash
# Quick development build
make dev

# Run comprehensive tests
make test

# Start documentation server with live reload
make livehtml

# Run sample simulation
cd src/supy/sample_run && suews-run .

# Convert legacy configurations
suews-convert to-yaml -i old_config -o new_config.yml
```

### With Claude Code

Claude Code provides intelligent assistance for:
- **Code Development**: Writing, refactoring, optimising SUEWS code
- **Debugging**: Identifying and fixing issues across Fortran/Python codebase
- **Testing**: Creating comprehensive test suites and interpreting results
- **Documentation**: Maintaining academic-standard documentation
- **Research**: Data analysis, visualisation, and publication preparation

## Academic Research Features

### Publication Workflow
- **Version Control**: All changes tracked in Git
- **Reproducibility**: Consistent environment across machines
- **Documentation**: Automatic generation of research documentation
- **British Standards**: All outputs follow British English conventions
- **Citation Management**: Integrated with academic citation tools

### Data Management
- **Test Data**: Extensive meteorological datasets included
- **Benchmark Data**: Reference results for validation
- **Output Organisation**: Structured directories for research outputs
- **Configuration Management**: Version-controlled SUEWS configurations

## Security

### Configuration Management
- **Development Settings**: Optional configuration in `.env` file
- **Build Configuration**: Parallel job settings and optimisation
- **Logging**: Configurable log levels for debugging
- **Access Control**: Read-only mounts for sensitive data

### Resource Protection
- **Memory Limits**: Prevent system resource exhaustion
- **Process Limits**: Control concurrent operations
- **Network Isolation**: Optional network restriction
- **File System**: Isolated container environment

## Troubleshooting

### Build Issues

```bash
# Check Docker installation
docker version

# Verify Claude Code Sandbox
claude-sandbox --version

# Check configuration syntax
cat claude-sandbox.config.json | jq .

# Clean and rebuild
docker system prune -f
./claude-dev/setup-claude-dev.sh
```

### Runtime Issues

```bash
# Check container logs
docker logs $(docker ps -q --filter ancestor=suews-research)

# Access container shell
docker exec -it $(docker ps -q --filter ancestor=suews-research) bash

# Check environment
micromamba list -n suews-dev

# Verify SUEWS installation
python -c "import supy; print(supy.__version__)"
```

### Performance Issues

```bash
# Monitor resource usage
docker stats

# Check memory usage inside container
free -h

# Profile application performance
python -m cProfile -o profile.stats your_script.py
```

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

*This Docker integration is optimised for academic research and follows British academic standards. For technical support, please use the GitHub issue tracker.*