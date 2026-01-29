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
# Download URLs
# Update these URLs after uploading to Zenodo/etc
#################################################

# Pre-trained model (~1.9 GB)
# Original source: TSignal authors via Dropbox
MODEL_URL="https://www.dropbox.com/s/lfuleg9470s7nqx/deployment_sep_pe_swa_extra_inpemb_on_gen_best_eval_only_dec.pth?dl=1"

# Singularity container (~557 MB)
# TODO: Update this URL after uploading to Zenodo
CONTAINER_URL="PLACEHOLDER_URL"

#################################################
# Download Pre-trained Model
#################################################

MODEL_DIR="model"
MODEL_NAME="deployment_sep_pe_swa_extra_inpemb_on_gen_best_eval_only_dec.pth"
MODEL_PATH="${MODEL_DIR}/${MODEL_NAME}"

mkdir -p "$MODEL_DIR"

if [[ -f "$MODEL_PATH" ]]; then
    echo "Model already exists: $MODEL_PATH"
    echo "  Size: $(du -h "$MODEL_PATH" | cut -f1)"
else
    echo "Downloading pre-trained model (~1.9 GB)..."
    echo "  This may take several minutes..."
    echo ""

    if [[ "$DOWNLOADER" == "wget" ]]; then
        wget --progress=bar:force -O "$MODEL_PATH" "$MODEL_URL"
    else
        curl -L --progress-bar -o "$MODEL_PATH" "$MODEL_URL"
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
# Download Singularity Container
#################################################

CONTAINER_DIR="container"
CONTAINER_PATH="${CONTAINER_DIR}/tsignal.sif"

mkdir -p "$CONTAINER_DIR"

if [[ -f "$CONTAINER_PATH" ]]; then
    echo "Container already exists: $CONTAINER_PATH"
    echo "  Size: $(du -h "$CONTAINER_PATH" | cut -f1)"
elif [[ "$CONTAINER_URL" == "PLACEHOLDER_URL" ]]; then
    echo "==========================================="
    echo "Container download URL not yet configured"
    echo "==========================================="
    echo ""
    echo "The container file needs to be downloaded manually."
    echo ""
    echo "Please contact the repository maintainer or check"
    echo "the GitHub releases page for the container file."
    echo ""
    echo "Once downloaded, place it at: $CONTAINER_PATH"
    echo ""
else
    echo "Downloading Singularity container (~557 MB)..."
    echo "  This may take a few minutes..."
    echo ""

    if [[ "$DOWNLOADER" == "wget" ]]; then
        wget --progress=bar:force -O "$CONTAINER_PATH" "$CONTAINER_URL"
    else
        curl -L --progress-bar -o "$CONTAINER_PATH" "$CONTAINER_URL"
    fi

    if [[ -f "$CONTAINER_PATH" ]]; then
        echo ""
        echo "Container downloaded successfully!"
        echo "  Size: $(du -h "$CONTAINER_PATH" | cut -f1)"
    else
        echo "ERROR: Container download failed!"
        exit 1
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
