#!/bin/bash
#
#SBATCH -p cpu
#SBATCH -n 1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=tsignal_%J_stdout.txt
#SBATCH --error=tsignal_%J_stderr.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=tsignal
#SBATCH --mail-user=YOUR_EMAIL@university.edu
#SBATCH --mail-type=ALL
#
#################################################
# TSignal Signal Peptide Prediction - SLURM Script
#
# Usage:
#   1. Edit INPUT_FASTA and OUTPUT_DIR below
#   2. Update email address above
#   3. Submit: sbatch tsignal_slurm.sbatch
#
# Author: Jorge L. Perez-Moreno, Ph.D. (jperezmoreno@umass.edu)
#################################################

#################################################
### CONFIGURE THESE VARIABLES
#################################################

# Input FASTA file (absolute path recommended)
INPUT_FASTA="sequences.fasta"

# Output directory
OUTPUT_DIR="output"

# SP probability threshold for removal (0.9 = high confidence)
THRESHOLD="0.9"

# Conda environment (if needed)
# CONDA_ENV="base"

#################################################
### SETUP
#################################################

echo "==========================================="
echo "TSignal Signal Peptide Prediction"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "==========================================="
echo ""

# Get script directory (where the toolkit is installed)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "$SCRIPT_DIR")"

echo "Toolkit directory: ${TOOLKIT_DIR}"
echo "Input: ${INPUT_FASTA}"
echo "Output: ${OUTPUT_DIR}"
echo "Threshold: ${THRESHOLD}"
echo ""

# Activate conda if needed
# source ~/.bashrc
# conda activate ${CONDA_ENV}

# Check Apptainer/Singularity
if command -v apptainer &> /dev/null; then
    RUNTIME="apptainer"
elif command -v singularity &> /dev/null; then
    RUNTIME="singularity"
else
    echo "ERROR: Neither apptainer nor singularity found"
    exit 1
fi
echo "Container runtime: ${RUNTIME}"

#################################################
### RUN TSIGNAL
#################################################

# Check input file
if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "ERROR: Input file not found: ${INPUT_FASTA}"
    exit 1
fi

echo ""
echo "Input sequences: $(grep -c '^>' ${INPUT_FASTA})"
echo ""

# Run batch processing
echo "Running TSignal batch processing..."
python3 "${SCRIPT_DIR}/tsignal_batch.py" \
    -i "${INPUT_FASTA}" \
    -o "${OUTPUT_DIR}" \
    --threshold "${THRESHOLD}" \
    --verbose

#################################################
### SUMMARY
#################################################

if [[ $? -eq 0 ]]; then
    echo ""
    echo "==========================================="
    echo "SUCCESS: TSignal processing complete"
    echo "==========================================="
    echo ""
    echo "Output files:"
    ls -lh "${OUTPUT_DIR}/"
    echo ""
    echo "Summary:"
    cat "${OUTPUT_DIR}/processing_summary.txt"
else
    echo ""
    echo "==========================================="
    echo "ERROR: TSignal processing failed"
    echo "==========================================="
    exit 1
fi

echo ""
echo "Finished: $(date)"
echo "==========================================="
