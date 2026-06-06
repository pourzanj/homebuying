#!/bin/bash
set -euo pipefail

# Docker-based local runner for daily_fit.R
# This replicates the GitHub Action using local Docker

# Load environment variables from .env file if it exists
if [ -f .env ]; then
    echo "Loading environment variables from .env file..."
    set -a
    . ./.env
    set +a
else
    echo "Warning: .env file not found. Using system environment variables."
fi

IMAGE_NAME="${HOMEBUYING_IMAGE:-homebuying:local}"

# Build the Docker image locally when it is missing, or when explicitly requested.
if [ "${REBUILD_IMAGE:-0}" = "1" ] || ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
    echo "Building local Docker image..."
    docker build -t "$IMAGE_NAME" .
else
    echo "Using existing Docker image $IMAGE_NAME. Set REBUILD_IMAGE=1 to rebuild."
fi

# Check if RIINGO_TOKEN is set (optional - only needed if using Tiingo endpoint)
if [ -z "${RIINGO_TOKEN:-}" ]; then
    echo "Note: RIINGO_TOKEN not set. Using Yahoo Finance endpoint (default)."
    echo "To use Tiingo instead, set RIINGO_TOKEN in .env file (see env.example)."
else
    echo "RIINGO_TOKEN found. Available for Tiingo endpoint if configured."
fi

# Set compute source if not already set
if [ -z "${COMPUTE_SOURCE:-}" ]; then
    export COMPUTE_SOURCE="docker-local"
fi

echo "Running daily_fit.R in Docker container (compute source: $COMPUTE_SOURCE)..."
docker run --rm \
    -e RIINGO_TOKEN="${RIINGO_TOKEN:-}" \
    -e COMPUTE_SOURCE="$COMPUTE_SOURCE" \
    "$IMAGE_NAME" \
    /bin/bash -euo pipefail -c 'tbb_lib=$(find /root/.cmdstan -path "*/stan/lib/stan_math/lib/tbb/libtbb.so.2" -print | sort -V | tail -n 1); if [ -z "$tbb_lib" ]; then echo "Could not find CmdStan TBB library under /root/.cmdstan" >&2; exit 1; fi; tbb_dir=$(dirname "$tbb_lib"); export LD_LIBRARY_PATH="$tbb_dir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"; R -e "renv::restore()"; Rscript scripts/daily_fit.R'

echo "Daily fit completed!"
