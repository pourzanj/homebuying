#!/bin/bash

# Local runner for daily_fit.R
# This script replicates what the GitHub Action does locally

echo "Setting up environment for daily_fit.R..."

# Load environment variables from .env file if it exists
if [ -f .env ]; then
    echo "Loading environment variables from .env file..."
    export $(grep -v '^#' .env | xargs)
else
    echo "Warning: .env file not found. Using system environment variables."
fi

# Check if R is installed
if ! command -v R &> /dev/null; then
    echo "Error: R is not installed. Please install R first."
    echo "Visit: https://cran.r-project.org/"
    exit 1
fi

# Check if RIINGO_TOKEN is set (optional - only needed if using Tiingo endpoint)
if [ -z "$RIINGO_TOKEN" ]; then
    echo "Note: RIINGO_TOKEN not set. Using Yahoo Finance endpoint (default)."
    echo "To use Tiingo instead, set RIINGO_TOKEN in .env file (see env.example)."
else
    echo "RIINGO_TOKEN found. Available for Tiingo endpoint if configured."
fi

# Set compute source if not already set in .env
if [ -z "$COMPUTE_SOURCE" ]; then
    export COMPUTE_SOURCE="local"
fi

echo "Compute source: $COMPUTE_SOURCE"

# Restore R dependencies
echo "Restoring R dependencies..."
R -e "renv::restore()"

# Install CmdStan if not already installed
echo "Installing CmdStan..."
R -e "cmdstanr::install_cmdstan()"

# Run the daily fit script
echo "Running daily_fit.R..."
Rscript scripts/daily_fit.R

echo "Daily fit completed!"
