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

    # Create start script
    cat > "$SUEWS_ROOT/start-claude-dev.sh" << 'EOF'
#!/bin/bash
# Quick start script for SUEWS Claude Code development

cd "$(dirname "$0")"

REBUILD_FLAG=false
OTHER_ARGS=()
BASE_CONTEXT=$(git rev-parse --abbrev-ref HEAD) # Default to current branch
CUSTOM_IMAGE="suews-claude-dev:latest"

# Separate --rebuild from other flags and determine the base context for the sandbox
for arg in "$@"; do
  if [[ "$arg" == "--rebuild" ]]; then
    REBUILD_FLAG=true
  else
    OTHER_ARGS+=("$arg")
  fi
done

# Find the user-specified context, if any
for i in "${!OTHER_ARGS[@]}"; do
  if [[ "${OTHER_ARGS[$i]}" == "--branch" || "${OTHER_ARGS[$i]}" == "-b" || "${OTHER_ARGS[$i]}" == "--remote-branch" ]]; then
    BASE_CONTEXT="${OTHER_ARGS[$i+1]}"
    break
  elif [[ "${OTHER_ARGS[$i]}" == "--pr" ]]; then
    BASE_CONTEXT="PR #${OTHER_ARGS[$i+1]}"
    break
  fi
done

# Build or rebuild the custom Docker image
build_custom_image() {
  echo "🏗️  Building custom SUEWS Claude Code Docker image..."
  docker build -f claude-dev/Dockerfile.claude-dev -t "$CUSTOM_IMAGE" . || {
    echo "❌ Failed to build Docker image"
    exit 1
  }
  echo "✅ Custom Docker image built: $CUSTOM_IMAGE"
}

# Check if custom image exists or rebuild is requested
if [[ "$REBUILD_FLAG" == true ]] || ! docker image inspect "$CUSTOM_IMAGE" >/dev/null 2>&1; then
  if [[ "$REBUILD_FLAG" == true ]]; then
    echo "🗑️  --rebuild flag detected. Forcing a full rebuild..."
    # Remove all possible images
    docker rmi -f claude-code-sandbox:latest > /dev/null 2>&1 || true
    docker rmi -f "$CUSTOM_IMAGE" > /dev/null 2>&1 || true
    docker rmi -f suews-claude-dev:test > /dev/null 2>&1 || true
    # Also clean up any claude-sandbox cached images
    claude-sandbox clean --images > /dev/null 2>&1 || true
    echo "✅ Images removed. Building fresh image..."
  else
    echo "🔍 Custom Docker image not found. Building it now..."
  fi
  build_custom_image
fi

echo "🚀 Starting SUEWS Claude Code development environment..."
echo "🌿 Base context: $BASE_CONTEXT"
echo "📁 Working directory: $(pwd)"
echo "🐳 Using custom image: $CUSTOM_IMAGE"
echo "🐳 Container config: ./claude-dev/claude-sandbox.config.json"
if [ ${#OTHER_ARGS[@]} -ne 0 ]; then
    echo " passing additional arguments: ${OTHER_ARGS[*]}"
fi
echo ""

export COPYFILE_DISABLE=1

# Start Claude Code Sandbox with custom image
CONFIG_FILE="./claude-dev/claude-sandbox.config.json"

# Check if jq is installed
if ! command -v jq > /dev/null; then
  echo "⚠️  jq is not installed. Mount paths with '~' may not work." >&2
  echo "   Attempting to start anyway..." >&2
  # Pass the custom image to claude-sandbox
  claude-sandbox start -c "$CONFIG_FILE" --image "$CUSTOM_IMAGE" "${OTHER_ARGS[@]}"
else
  # Use jq to dynamically replace ~ with the user's home directory
  # This ensures that mounts for .ssh and .gitconfig work correctly
  # Also inject the custom image into the config
  jq ". + {\"image\": \"$CUSTOM_IMAGE\"} | (.mounts[] | select(.source | startswith(\"~\"))).source |= \"$HOME\" + (. | ltrimstr(\"~\"))" "$CONFIG_FILE" | \
  claude-sandbox start -c - "${OTHER_ARGS[@]}"
fi
EOF
    chmod +x "$SUEWS_ROOT/start-claude-dev.sh"

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
        "$SUEWS_ROOT/start-claude-dev.sh"
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
    echo "1. Run: ./start-claude-dev.sh"
    echo "2. Read: claude-dev/README.md for detailed usage"
    echo ""
    log_info "Happy coding with Claude! 🤖"
}

# Run main function
main "$@"