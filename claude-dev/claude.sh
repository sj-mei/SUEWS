#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.

# --- Configuration ---
# Default parent directory for all workspaces. Can be overridden with the LOCATION environment variable.
LOCATION="${LOCATION:-$HOME/claude-suews-workspace}"
# The root of the SUEWS project, determined by the script's location.
SUEWS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# List of directories/files to exclude when creating a new workspace.
RSYNC_EXCLUDES=(
  --exclude=".github" --exclude=".DS_Store" --exclude="build"
  --exclude="dist" --exclude="*.egg-info" --exclude="wheelhouse"
  --exclude="**/__pycache__" --exclude="*.pyc" --exclude="claude-dev/workspace"
  --exclude="*.log" --exclude="docs/build" --exclude="Release"
  --exclude=".vscode" --exclude=".idea" --exclude=".venv"
)

# --- Helper Functions ---
print_usage() {
  cat <<EOF
SUEWS Claude Development Environment Manager

This script helps create and manage isolated, parallel development workspaces.

Usage:
  ./claude.sh <command> [arguments]

Commands:
  dev [name]            Create a new development workspace.
                        - If [name] is provided, creates a named workspace (e.g., SUEWS-my-feature).
                        - If no name is provided, creates a timestamped workspace.

  start [name] [args]   Start the Docker container for a given workspace.
                        - If [name] is omitted, an interactive selection menu is shown.
                        - Any additional [args] are passed to the start-claude-dev.sh script.

  stop <name>           Stop the Docker container for a given workspace.

  list                  List all available SUEWS workspaces in:
                        ${LOCATION}

  clean <name> [force]  Clean (delete) a specific workspace.
                        - Prompts for confirmation unless 'force' is provided as the second argument.

  clean-all [force]     Clean (delete) ALL workspaces.
                        - DANGEROUS: This removes the entire workspace parent directory:
                          ${LOCATION}

  help                  Show this help message.

Examples:
  ./claude-dev/claude.sh dev feature-x
  ./claude-dev/claude.sh start feature-x
  ./claude-dev/claude.sh list
  ./claude-dev/claude.sh stop feature-x
  ./claude-dev/claude.sh clean feature-x

EOF
}

# --- Command Implementations ---

cmd_dev() {
  local timestamp
  timestamp=$(date +%Y%m%d-%H%M%S)
  local copy_name="${1:-$timestamp}"
  local target_dir="${LOCATION}/SUEWS-${copy_name}"

  echo "🚀 Creating a parallel SUEWS workspace with git history..."
  echo "   Source:      ${SUEWS_ROOT}"
  echo "   Destination: ${target_dir}"
  echo ""

  mkdir -p "${LOCATION}"

  if [ -d "${target_dir}" ]; then
    echo "❌ Error: Target directory '${target_dir}' already exists." >&2
    echo "   Please choose a different name." >&2
    exit 1
  fi

  local remote_url
  remote_url=$(cd "${SUEWS_ROOT}" && git remote get-url origin)
  if [ -z "$remote_url" ]; then
    echo "❌ Error: Could not determine git remote URL for 'origin'." >&2
    exit 1
  fi

  local source_branch
  source_branch=$(cd "${SUEWS_ROOT}" && git rev-parse --abbrev-ref HEAD)
  if [ -z "$source_branch" ]; then
    echo "❌ Error: Could not determine current git branch." >&2
    exit 1
  fi

  echo "📦 Cloning repository from '${remote_url}' (branch: ${source_branch})..."
  echo "   This will create a clean copy with committed changes only."
  git clone --recurse-submodules --branch "${source_branch}" "${remote_url}" "${target_dir}"

  echo ""
  echo "🔧 Running workspace setup to generate launch scripts..."
  cd "${target_dir}"
  chmod +x claude-dev/setup-claude-dev.sh
  ./claude-dev/setup-claude-dev.sh

  echo ""
  echo "✅ Workspace with git history created successfully!"
  echo "📂 Location: ${target_dir}"
  echo ""
  echo "💡 This workspace is a clean clone of the '${source_branch}' branch."
  echo "🚀 To start working safely, you can create a new feature branch:"
  echo "   cd \"${target_dir}\" && git checkout -b feature/${copy_name}"
}

cmd_start() {
  local name="$1"
  # All other arguments are passed through to the underlying script
  # The first shift removes the 'name' if it was provided
  [ -n "$1" ] && shift
  local extra_args=("$@")

  if [ -z "$name" ]; then
    echo "▶️  No workspace specified. Please choose one to start:"

    local workspaces_raw
    workspaces_raw=($(ls -d "${LOCATION}/SUEWS-"* 2>/dev/null))

    if [ ${#workspaces_raw[@]} -eq 0 ]; then
      echo "   No workspaces found."
      exit 0
    fi

    local workspaces=()
    for ws in "${workspaces_raw[@]}"; do
      workspaces+=("$(basename "$ws" | sed 's/SUEWS-//')")
    done

    PS3="Enter number: "
    select choice in "${workspaces[@]}"; do
      if [[ -n "$choice" ]]; then
        name="$choice"
        break
      else
        echo "Invalid selection. Please try again." >&2
      fi
    done
    echo "" # Add a newline for better formatting
  fi

  local target_dir="${LOCATION}/SUEWS-${name}"

  if [ ! -d "${target_dir}" ]; then
    echo "❌ Error: Workspace 'SUEWS-${name}' not found in '${LOCATION}'" >&2
    exit 1
  fi

  echo "▶️  Starting Claude Code sandbox for 'SUEWS-${name}'..."
  cd "${target_dir}"
  # Pass any remaining arguments to the script
  ./start-claude-dev.sh "${extra_args[@]}"
}

cmd_stop() {
  local name="$1"
  local target_dir="${LOCATION}/SUEWS-${name}"

  if [ -z "$name" ]; then
    echo "❌ Error: Please specify which workspace to stop." >&2
    echo "   Usage: ./claude.sh stop <name>" >&2
    exit 1
  fi

  if [ ! -d "${target_dir}" ]; then
    echo "❌ Error: Workspace 'SUEWS-${name}' not found in '${LOCATION}'" >&2
    exit 1
  fi

  echo "⏹️  Stopping Claude Code sandbox for 'SUEWS-${name}'..."
  cd "${target_dir}"
  ./stop-claude-dev.sh
}

cmd_list() {
  echo "📋 Available SUEWS workspaces in: ${LOCATION}"
  echo ""
  if ! ls -d "${LOCATION}/SUEWS-"* 2>/dev/null | sed "s|.*/SUEWS-||g"; then
    echo "   No workspaces found."
  fi
  echo ""
}

cmd_clean() {
  local name="$1"
  local force="$2"
  local target_dir="${LOCATION}/SUEWS-${name}"

  if [ -z "$name" ]; then
    echo "❌ Error: Please specify which workspace to clean." >&2
    echo "   Usage: ./claude.sh clean <name>" >&2
    exit 1
  fi

  if [ ! -d "${target_dir}" ]; then
    echo "ℹ️  Workspace 'SUEWS-${name}' not found in '${LOCATION}'. Nothing to do."
    exit 0
  fi

  echo "🧹 Preparing to clean workspace 'SUEWS-${name}'..."
  if [ "$force" == "force" ]; then
    echo "📂 Removing workspace (forced): ${target_dir}"
    rm -rf "${target_dir}"
    echo "✅ Workspace removed."
  else
    read -rp "❓ Are you sure you want to delete '${target_dir}'? This action cannot be undone. [y/N] " REPLY
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      echo "📂 Removing workspace: ${target_dir}"
      rm -rf "${target_dir}"
      echo "✅ Workspace removed."
    else
      echo "🛑 Aborted by user."
    fi
  fi
}

cmd_clean_all() {
  local force="$1"

  if [ ! -d "${LOCATION}" ] || [ -z "$(ls -A "${LOCATION}")" ]; then
    echo "ℹ️  Workspace directory '${LOCATION}' is empty or does not exist. Nothing to do."
    exit 0
  fi

  echo "🧹 Preparing to clean ALL workspaces in: ${LOCATION}"
  echo "The following workspaces will be PERMANENTLY DELETED:"
  ls -d "${LOCATION}/SUEWS-"* 2>/dev/null | sed "s|${LOCATION}/||g" | sed 's/^/  - /' || echo "   (No workspaces found to delete)"
  echo ""

  if [ "$force" == "force" ]; then
    echo "📂 Removing all workspaces (forced)..."
    rm -rf "${LOCATION}"
    echo "✅ All workspaces removed."
  else
    read -rp "❓ Are you sure you want to delete ALL workspaces? This action cannot be undone. [y/N] " REPLY
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      echo "📂 Removing all workspaces..."
      rm -rf "${LOCATION}"
      echo "✅ All workspaces removed."
    else
      echo "🛑 Aborted by user."
    fi
  fi
}

# --- Main Dispatcher ---
main() {
  local cmd="$1"
  if [ -z "$cmd" ]; then
    print_usage
    exit 1
  fi
  shift # remove command from arguments

  case "$cmd" in
    dev)
      cmd_dev "$@"
      ;;
    start)
      cmd_start "$@"
      ;;
    stop)
      cmd_stop "$@"
      ;;
    list)
      cmd_list "$@"
      ;;
    clean)
      cmd_clean "$@"
      ;;
    clean-all)
      cmd_clean_all "$@"
      ;;
    help|-h|--help)
      print_usage
      ;;
    *)
      echo "❌ Error: Unknown command '$cmd'" >&2
      echo ""
      print_usage
      exit 1
      ;;
  esac
}

main "$@"