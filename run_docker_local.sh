#!/bin/bash

# Docker-based local runner for daily_fit.R
# This replicates the GitHub Action using local Docker

echo "Building local Docker image..."

# Build the Docker image locally
docker build -t homebuying:local .

# Load environment variables from .env file if it exists
if [ -f .env ]; then
    echo "Loading environment variables from .env file..."
    export $(grep -v '^#' .env | xargs)
else
    echo "Warning: .env file not found. Using system environment variables."
fi

# Check if RIINGO_TOKEN is set (optional - only needed if using Tiingo endpoint)
if [ -z "$RIINGO_TOKEN" ]; then
    echo "Note: RIINGO_TOKEN not set. Using Yahoo Finance endpoint (default)."
    echo "To use Tiingo instead, set RIINGO_TOKEN in .env file (see env.example)."
else
    echo "RIINGO_TOKEN found. Available for Tiingo endpoint if configured."
fi

# Set compute source if not already set
if [ -z "$COMPUTE_SOURCE" ]; then
    export COMPUTE_SOURCE="docker-local"
fi

echo "Running daily_fit.R in Docker container (compute source: $COMPUTE_SOURCE)..."
docker run --rm \
    -e RIINGO_TOKEN="${RIINGO_TOKEN:-}" \
    -e COMPUTE_SOURCE=$COMPUTE_SOURCE \
    homebuying:local \
    /bin/bash -c "export LD_LIBRARY_PATH=/root/.cmdstan/cmdstan-2.38.0/stan/lib/stan_math/lib/tbb:\$LD_LIBRARY_PATH; R -e 'renv::restore()'; Rscript scripts/daily_fit.R"

echo "Daily fit completed!"
