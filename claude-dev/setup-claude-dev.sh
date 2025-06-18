#!/bin/bash
# SUEWS Claude Code Development Environment Setup
# 
# This script sets up a complete development environment for SUEWS
# with Claude Code integration for enhanced productivity.
#
# Usage: ./setup-claude-dev.sh
# 
# Prerequisites:
# - Docker or Podman
# - Node.js and npm
# - Git (for cloning dependencies)

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUEWS_ROOT="$(dirname "$SCRIPT_DIR")"
DOCKER_DIR="$SUEWS_ROOT/claude-dev"

# Colours for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

log_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

log_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

log_error() {
    echo -e "${RED}❌ $1${NC}"
}

log_header() {
    echo -e "${CYAN}🚀 $1${NC}"
}

# Check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check prerequisites
check_prerequisites() {
    log_header "Checking prerequisites..."
    
    local missing_deps=()
    
    # Check Docker or Podman
    if ! command_exists docker && ! command_exists podman; then
        missing_deps+=("docker or podman")
    fi
    
    # Check Node.js
    if ! command_exists node; then
        missing_deps+=("node.js")
    fi
    
    # Check npm
    if ! command_exists npm; then
        missing_deps+=("npm")
    fi
    
    # Check git
    if ! command_exists git; then
        missing_deps+=("git")
    fi
    
    if [ ${#missing_deps[@]} -ne 0 ]; then
        log_error "Missing dependencies: ${missing_deps[*]}"
        echo ""
        echo "Please install the missing dependencies:"
        echo "  macOS: brew install docker node git"
        echo "  Ubuntu: sudo apt-get install docker.io nodejs npm git"
        echo "  Fedora: sudo dnf install docker nodejs npm git"
        exit 1
    fi
    
    log_success "All prerequisites satisfied"
}

# Install Claude Code Sandbox
install_claude_sandbox() {
    log_header "Installing Claude Code Sandbox..."
    
    if command_exists claude-sandbox; then
        log_info "Claude Code Sandbox already installed"
        return
    fi
    
    log_info "Installing Claude Code and Sandbox tools..."
    npm install -g @anthropic-ai/claude-code @textcortex/claude-code-sandbox
    
    if command_exists claude-sandbox; then
        log_success "Claude Code Sandbox installed successfully"
    else
        log_error "Failed to install Claude Code Sandbox"
        exit 1
    fi
}

# Create directory structure
create_directories() {
    log_header "Creating directory structure..."
    
    # Create optional directories for mounts
    mkdir -p "$SUEWS_ROOT/notebooks"
    mkdir -p "$SUEWS_ROOT/configs" 
    mkdir -p "$SUEWS_ROOT/outputs"
    mkdir -p "$SUEWS_ROOT/external-data"
    
    log_success "Directory structure created"
}

# Create environment file template
create_env_template() {
    log_header "Creating environment configuration..."
    
    local env_file="$SUEWS_ROOT/.env"
    
    if [ -f "$env_file" ]; then
        log_info "Environment file already exists: $env_file"
        return
    fi
    
    cat > "$env_file" << 'EOF'
# SUEWS Development Environment Configuration
# 
# Optional configuration file for development settings

# === Development Configuration ===
# Set to 'true' to enable development mode features
SUEWS_DEV_MODE=true

# Number of parallel build jobs (adjust based on your system)
SUEWS_BUILD_JOBS=4

# === External Data Sources (Optional) ===
# ERA5 CDS API credentials (if using external ERA5 data)
# CDSAPI_URL=https://cds.climate.copernicus.eu/api/v2
# CDSAPI_KEY=your-cds-api-key

# === Logging Configuration ===
# Log level for development (DEBUG, INFO, WARNING, ERROR)
SUEWS_LOG_LEVEL=INFO
EOF
    
    log_success "Environment template created: $env_file"
}

# Create convenience scripts
create_scripts() {
    log_header "Creating convenience scripts..."
    
    # Create run script
    cat > "$SUEWS_ROOT/run-claude-dev.sh" << 'EOF'
#!/bin/bash
# Quick start script for SUEWS Claude Code development

cd "$(dirname "$0")"

echo "🚀 Starting SUEWS Claude Code development environment..."
echo "📁 Working directory: $(pwd)"
echo "🐳 Container config: ./claude-dev/claude-sandbox.config.json"
echo ""

# Start Claude Code Sandbox
claude-sandbox start -c ./claude-dev/claude-sandbox.config.json
EOF
    chmod +x "$SUEWS_ROOT/run-claude-dev.sh"
    
    # Create stop script
    cat > "$SUEWS_ROOT/stop-claude-dev.sh" << 'EOF'
#!/bin/bash
# Stop SUEWS Claude Code development environment

echo "🛑 Stopping SUEWS Claude Code development environment..."
claude-sandbox stop
echo "✅ Environment stopped"
EOF
    chmod +x "$SUEWS_ROOT/stop-claude-dev.sh"
    
    # Create cleanup script
    cat > "$SUEWS_ROOT/cleanup-claude-dev.sh" << 'EOF'
#!/bin/bash
# Clean up SUEWS Claude Code development environment

echo "🧹 Cleaning up SUEWS Claude Code development environment..."
claude-sandbox clean --force
echo "✅ Cleanup complete"
EOF
    chmod +x "$SUEWS_ROOT/cleanup-claude-dev.sh"
    
    log_success "Convenience scripts created"
}


# Verify installation
verify_installation() {
    log_header "Verifying installation..."
    
    # Check files exist
    local required_files=(
        "$DOCKER_DIR/Dockerfile.claude-dev"
        "$DOCKER_DIR/claude-sandbox.config.json"
        "$SUEWS_ROOT/.env"
        "$SUEWS_ROOT/run-claude-dev.sh"
    )
    
    for file in "${required_files[@]}"; do
        if [ -f "$file" ]; then
            log_success "✓ $(basename "$file")"
        else
            log_error "✗ $(basename "$file") missing"
            return 1
        fi
    done
    
    # Check Claude Code Sandbox
    if command_exists claude-sandbox; then
        log_success "✓ claude-sandbox command available"
    else
        log_error "✗ claude-sandbox command not found"
        return 1
    fi
    
    log_success "Installation verification complete"
}

# Main setup function
main() {
    echo ""
    log_header "SUEWS Claude Code Development Environment Setup"
    echo ""
    
    check_prerequisites
    install_claude_sandbox
    create_directories
    create_env_template
    create_scripts
    verify_installation
    
    echo ""
    log_success "🎉 Setup complete!"
    echo ""
    echo "Next steps:"
    echo "1. Run: ./run-claude-dev.sh"
    echo "2. Read: claude-dev/README.md for detailed usage"
    echo ""
    log_info "Happy coding with Claude! 🤖"
}

# Run main function
main "$@"