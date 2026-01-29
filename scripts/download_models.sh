#!/bin/bash
#################################################
# Download TSignal Model and Container
#
# This script downloads the large files required
# to run TSignal that are not included in the
# GitHub repository due to size constraints.
#
# Author: Jorge L. Perez-Moreno, Ph.D.
#################################################

set -e

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$TOOLKIT_DIR"

echo "==========================================="
echo "TSignal Toolkit - Download Required Files"
echo "==========================================="
echo ""

# Check for wget or curl
if command -v wget &> /dev/null; then
    DOWNLOADER="wget"
elif command -v curl &> /dev/null; then
    DOWNLOADER="curl"
else
    echo "ERROR: Neither wget nor curl found. Please install one of them."
    exit 1
fi

echo "Using: $DOWNLOADER"
echo ""

#################################################
# Download Pre-trained Model
#################################################

MODEL_DIR="model"
MODEL_NAME="deployment_sep_pe_swa_extra_inpemb_on_gen_best_eval_only_dec.pth"
MODEL_PATH="${MODEL_DIR}/${MODEL_NAME}"
MODEL_URL="https://www.dropbox.com/s/lfuleg9470s7nqx/${MODEL_NAME}?dl=1"

mkdir -p "$MODEL_DIR"

if [[ -f "$MODEL_PATH" ]]; then
    echo "Model already exists: $MODEL_PATH"
    echo "  Size: $(du -h "$MODEL_PATH" | cut -f1)"
else
    echo "Downloading pre-trained model (~1.9 GB)..."
    echo "  URL: $MODEL_URL"
    echo "  Destination: $MODEL_PATH"
    echo ""

    if [[ "$DOWNLOADER" == "wget" ]]; then
        wget -O "$MODEL_PATH" "$MODEL_URL"
    else
        curl -L -o "$MODEL_PATH" "$MODEL_URL"
    fi

    if [[ -f "$MODEL_PATH" ]]; then
        echo ""
        echo "Model downloaded successfully!"
        echo "  Size: $(du -h "$MODEL_PATH" | cut -f1)"
    else
        echo "ERROR: Model download failed!"
        exit 1
    fi
fi

echo ""

#################################################
# Build Singularity Container
#################################################

CONTAINER_DIR="container"
CONTAINER_PATH="${CONTAINER_DIR}/tsignal.sif"
DEF_PATH="${CONTAINER_DIR}/tsignal.def"

mkdir -p "$CONTAINER_DIR"

if [[ -f "$CONTAINER_PATH" ]]; then
    echo "Container already exists: $CONTAINER_PATH"
    echo "  Size: $(du -h "$CONTAINER_PATH" | cut -f1)"
else
    echo "Building Singularity container (~557 MB)..."
    echo ""
    echo "NOTE: This requires either:"
    echo "  1. sudo privileges (for local build), or"
    echo "  2. A Sylabs account (for remote build)"
    echo ""

    # Check for singularity/apptainer
    if command -v apptainer &> /dev/null; then
        RUNTIME="apptainer"
    elif command -v singularity &> /dev/null; then
        RUNTIME="singularity"
    else
        echo "WARNING: Neither apptainer nor singularity found."
        echo ""
        echo "To install on Ubuntu/Debian:"
        echo "  sudo apt install apptainer"
        echo ""
        echo "Or on HPC systems:"
        echo "  module load singularity"
        echo ""
        echo "After installing, re-run this script or build manually:"
        echo "  sudo singularity build $CONTAINER_PATH $DEF_PATH"
        echo ""
        echo "Skipping container build for now..."
        RUNTIME=""
    fi

    if [[ -n "$RUNTIME" ]]; then
        echo "Found: $RUNTIME"
        echo ""

        # Try to build
        if [[ -f "$DEF_PATH" ]]; then
            echo "Attempting to build container..."
            echo "  Definition: $DEF_PATH"
            echo "  Output: $CONTAINER_PATH"
            echo ""

            # Try remote build first (doesn't require sudo)
            if $RUNTIME build --remote "$CONTAINER_PATH" "$DEF_PATH" 2>/dev/null; then
                echo ""
                echo "Container built successfully (remote)!"
            elif sudo $RUNTIME build "$CONTAINER_PATH" "$DEF_PATH" 2>/dev/null; then
                echo ""
                echo "Container built successfully (local with sudo)!"
            else
                echo ""
                echo "WARNING: Container build failed."
                echo ""
                echo "Please build manually using one of these methods:"
                echo ""
                echo "  Option 1 - Local build (requires sudo):"
                echo "    sudo $RUNTIME build $CONTAINER_PATH $DEF_PATH"
                echo ""
                echo "  Option 2 - Remote build (requires Sylabs account):"
                echo "    $RUNTIME build --remote $CONTAINER_PATH $DEF_PATH"
                echo ""
            fi
        else
            echo "ERROR: Definition file not found: $DEF_PATH"
        fi
    fi
fi

echo ""

#################################################
# Verify Installation
#################################################

echo "==========================================="
echo "Verifying Installation"
echo "==========================================="
echo ""

ERRORS=0

if [[ -f "$MODEL_PATH" ]]; then
    SIZE=$(du -h "$MODEL_PATH" | cut -f1)
    echo "✓ Model found: $MODEL_PATH ($SIZE)"
else
    echo "✗ Model NOT found: $MODEL_PATH"
    ERRORS=$((ERRORS + 1))
fi

if [[ -f "$CONTAINER_PATH" ]]; then
    SIZE=$(du -h "$CONTAINER_PATH" | cut -f1)
    echo "✓ Container found: $CONTAINER_PATH ($SIZE)"
else
    echo "✗ Container NOT found: $CONTAINER_PATH"
    echo "  -> Build with: sudo singularity build $CONTAINER_PATH $DEF_PATH"
    ERRORS=$((ERRORS + 1))
fi

echo ""

if [[ $ERRORS -eq 0 ]]; then
    echo "==========================================="
    echo "Setup Complete!"
    echo "==========================================="
    echo ""
    echo "Test with:"
    echo "  ./scripts/tsignal_batch.py -i examples/test_sequences.fasta -o test_output/"
    echo ""
else
    echo "==========================================="
    echo "Setup Incomplete - $ERRORS file(s) missing"
    echo "==========================================="
    echo ""
    echo "See instructions above to complete setup."
    exit 1
fi
