#!/bin/bash
#################################################
# TSignal Toolkit Setup Script
#
# This script copies/links required files from the
# TSignal repository to set up the toolkit.
#
# Usage: ./setup.sh [path_to_TSignal_dir]
#
# Author: Jorge L. Perez-Moreno, Ph.D. (jperezmoreno@umass.edu)
#################################################

set -e

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Default TSignal source directory
TSIGNAL_DIR="${1:-../TSignal}"

echo "==========================================="
echo "TSignal Toolkit Setup"
echo "==========================================="
echo ""
echo "Toolkit directory: $SCRIPT_DIR"
echo "TSignal source: $TSIGNAL_DIR"
echo ""

# Check TSignal directory exists
if [[ ! -d "$TSIGNAL_DIR" ]]; then
    echo "ERROR: TSignal directory not found: $TSIGNAL_DIR"
    echo ""
    echo "Usage: $0 [path_to_TSignal_dir]"
    echo ""
    echo "The TSignal directory should contain:"
    echo "  - tsignal.sif (Singularity container)"
    echo "  - main.py"
    echo "  - sp_data/ directory with model file"
    exit 1
fi

# Create directories
echo "Creating directory structure..."
mkdir -p container model source/sp_data output

# Copy container
echo ""
echo "Copying Singularity container..."
if [[ -f "$TSIGNAL_DIR/tsignal.sif" ]]; then
    cp -v "$TSIGNAL_DIR/tsignal.sif" container/
else
    echo "WARNING: tsignal.sif not found - you'll need to copy it manually"
fi

# Copy model
echo ""
echo "Copying pre-trained model..."
MODEL_NAME="deployment_sep_pe_swa_extra_inpemb_on_gen_best_eval_only_dec.pth"
if [[ -f "$TSIGNAL_DIR/sp_data/$MODEL_NAME" ]]; then
    cp -v "$TSIGNAL_DIR/sp_data/$MODEL_NAME" model/
else
    echo "WARNING: Model file not found - you'll need to download it"
    echo "  Download from: https://www.dropbox.com/s/lfuleg9470s7nqx/$MODEL_NAME"
fi

# Copy source files
echo ""
echo "Copying source files..."
for item in main.py sp6_dicts.bin tokenizer_patch.py; do
    if [[ -f "$TSIGNAL_DIR/$item" ]]; then
        cp -v "$TSIGNAL_DIR/$item" source/
    fi
done

for dir in models utils misc; do
    if [[ -d "$TSIGNAL_DIR/$dir" ]]; then
        cp -rv "$TSIGNAL_DIR/$dir" source/
    fi
done

# Copy sp_data auxiliary files (not the large model)
echo ""
echo "Copying sp_data auxiliary files..."
for item in bert_tuning.py bert_tuning_tnmt.py data_utils.py create_test_files.py __init__.py sp6_data; do
    if [[ -e "$TSIGNAL_DIR/sp_data/$item" ]]; then
        cp -rv "$TSIGNAL_DIR/sp_data/$item" source/sp_data/
    fi
done

# Create symlink to model in sp_data
echo ""
echo "Creating model symlink..."
ln -sf "../../model/$MODEL_NAME" source/sp_data/

# Make scripts executable
echo ""
echo "Making scripts executable..."
chmod +x scripts/*.py scripts/*.sh

# Verify setup
echo ""
echo "==========================================="
echo "Verifying installation..."
echo "==========================================="
python3 scripts/tsignal_batch.py --check

echo ""
echo "==========================================="
echo "Setup complete!"
echo "==========================================="
echo ""
echo "To test the installation:"
echo "  ./scripts/tsignal_batch.py -i examples/test_sequences.fasta -o output/"
echo ""
echo "For more information, see README.md"
