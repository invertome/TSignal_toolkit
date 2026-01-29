#!/bin/bash
#################################################
# TSignal Prediction Script
# Run TSignal signal peptide prediction only
#
# Usage: ./tsignal_predict.sh input.fasta output.csv
#
# Author: Jorge L. Perez-Moreno, Ph.D. (jperezmoreno@umass.edu)
#################################################

set -e

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "$SCRIPT_DIR")"

# Configuration
CONTAINER="${TOOLKIT_DIR}/container/tsignal.sif"
SOURCE_DIR="${TOOLKIT_DIR}/source"
MODEL_NAME="deployment_sep_pe_swa_extra_inpemb_on_gen_best_eval_only_dec.pth"

# Check arguments
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <input.fasta> <output.csv> [--verbose]"
    echo ""
    echo "Arguments:"
    echo "  input.fasta   Input FASTA file with protein sequences"
    echo "  output.csv    Output CSV file for predictions"
    echo "  --verbose     Print detailed output"
    exit 1
fi

INPUT_FASTA="$1"
OUTPUT_CSV="$2"
VERBOSE=""
[[ "$3" == "--verbose" ]] && VERBOSE="--verbouse"

# Validate input
if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "ERROR: Input file not found: $INPUT_FASTA"
    exit 1
fi

# Check container
if [[ ! -f "$CONTAINER" ]]; then
    echo "ERROR: Container not found: $CONTAINER"
    echo "  -> Copy from: ../TSignal/tsignal.sif"
    exit 1
fi

# Check source directory
if [[ ! -d "$SOURCE_DIR" ]]; then
    echo "ERROR: Source directory not found: $SOURCE_DIR"
    echo "  -> Copy TSignal source files to source/"
    exit 1
fi

# Detect container runtime
if command -v apptainer &> /dev/null; then
    RUNTIME="apptainer"
elif command -v singularity &> /dev/null; then
    RUNTIME="singularity"
else
    echo "ERROR: Neither 'apptainer' nor 'singularity' found"
    exit 1
fi

echo "==========================================="
echo "TSignal Signal Peptide Prediction"
echo "==========================================="
echo "Input: $INPUT_FASTA"
echo "Output: $OUTPUT_CSV"
echo "Container: $RUNTIME"
echo "Started: $(date)"
echo ""

# Prepare working directory
WORK_DIR="$SOURCE_DIR/sp_data"
mkdir -p "$WORK_DIR"

# Copy input file
INPUT_NAME=$(basename "$INPUT_FASTA")
cp "$INPUT_FASTA" "$WORK_DIR/$INPUT_NAME"

# Link model if needed
if [[ ! -f "$WORK_DIR/$MODEL_NAME" ]]; then
    MODEL_PATH="${TOOLKIT_DIR}/model/${MODEL_NAME}"
    if [[ -f "$MODEL_PATH" ]]; then
        ln -sf "$MODEL_PATH" "$WORK_DIR/$MODEL_NAME"
    else
        echo "ERROR: Model not found: $MODEL_PATH"
        exit 1
    fi
fi

# Run TSignal
OUTPUT_NAME=$(basename "$OUTPUT_CSV")

echo "Running TSignal..."
$RUNTIME run \
    --bind "$SOURCE_DIR:/app/TSignal" \
    "$CONTAINER" \
    --test_seqs "sp_data/$INPUT_NAME" \
    --test_mdl "$MODEL_NAME" \
    --tune_bert \
    --train_only_decoder \
    --output_file "sp_data/$OUTPUT_NAME" \
    $VERBOSE

# Move output
if [[ -f "$WORK_DIR/$OUTPUT_NAME" ]]; then
    mv "$WORK_DIR/$OUTPUT_NAME" "$OUTPUT_CSV"
    echo ""
    echo "Predictions saved to: $OUTPUT_CSV"
else
    echo "ERROR: TSignal did not produce output file"
    exit 1
fi

# Cleanup
rm -f "$WORK_DIR/$INPUT_NAME"

echo ""
echo "Finished: $(date)"
echo "==========================================="
