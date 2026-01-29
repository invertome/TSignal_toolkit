#!/usr/bin/env python3
"""
TSignal Batch Processing Script
================================
All-in-one script for signal peptide prediction and removal.

Usage:
    ./tsignal_batch.py -i sequences.fasta -o output/
    ./tsignal_batch.py -i sequences.fasta -o output/ --threshold 0.8
    ./tsignal_batch.py --check  # Verify installation

Author: Jorge L. Perez-Moreno, Ph.D. (jperezmoreno@umass.edu)
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from datetime import datetime


# ============================================================================
# Configuration
# ============================================================================

TOOLKIT_DIR = Path(__file__).resolve().parent.parent
CONTAINER_PATH = TOOLKIT_DIR / "container" / "tsignal.sif"
MODEL_NAME = "deployment_sep_pe_swa_extra_inpemb_on_gen_best_eval_only_dec.pth"
MODEL_PATH = TOOLKIT_DIR / "model" / MODEL_NAME
SOURCE_DIR = TOOLKIT_DIR / "source"

# Signal peptide type mapping
SP_TYPES = {
    'S': 'Sec/SPase_I',    # Standard secretory
    'T': 'Tat/SPase_I',    # Twin-arginine translocation
    'L': 'Sec/SPase_II',   # Lipoprotein
    'W': 'Tat/SPase_II',   # Tat lipoprotein
    'P': 'Sec/SPase_IV',   # Prepilin-like
    'I': 'No_SP',          # No signal peptide
    'O': 'Other',          # Post-cleavage
    'M': 'Membrane',       # Transmembrane
    'E': 'End'             # Sequence end
}


# ============================================================================
# Utility Functions
# ============================================================================

def log(msg, level="INFO"):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {level}: {msg}", file=sys.stderr)


def check_installation():
    """Verify all required files are present."""
    errors = []

    # Check container
    if not CONTAINER_PATH.exists():
        errors.append(f"Container not found: {CONTAINER_PATH}")
        errors.append("  -> Copy from: ../TSignal/tsignal.sif")
    else:
        log(f"Container found: {CONTAINER_PATH} ({CONTAINER_PATH.stat().st_size / 1e6:.1f} MB)")

    # Check model
    if not MODEL_PATH.exists():
        errors.append(f"Model not found: {MODEL_PATH}")
        errors.append("  -> Copy from: ../TSignal/sp_data/" + MODEL_NAME)
    else:
        log(f"Model found: {MODEL_PATH} ({MODEL_PATH.stat().st_size / 1e9:.2f} GB)")

    # Check source files
    required_source = ["main.py", "models", "utils"]
    for item in required_source:
        path = SOURCE_DIR / item
        if not path.exists():
            errors.append(f"Source not found: {path}")

    if not errors:
        log("All required files present!")

        # Check singularity/apptainer
        for cmd in ["apptainer", "singularity"]:
            try:
                result = subprocess.run([cmd, "--version"], capture_output=True, text=True)
                if result.returncode == 0:
                    log(f"Container runtime: {cmd} {result.stdout.strip()}")
                    break
            except FileNotFoundError:
                continue
        else:
            errors.append("Neither 'apptainer' nor 'singularity' found in PATH")

    if errors:
        log("Installation check FAILED", "ERROR")
        for err in errors:
            print(f"  {err}", file=sys.stderr)
        return False

    log("Installation check PASSED")
    return True


def read_fasta(fasta_path):
    """Read FASTA file and return dict of {header: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]  # Remove '>'
                current_seq = []
            elif line:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def write_fasta(sequences, output_path, line_width=60):
    """Write sequences dict to FASTA file."""
    with open(output_path, 'w') as f:
        for header, seq in sequences.items():
            f.write(f'>{header}\n')
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')


# ============================================================================
# TSignal Prediction
# ============================================================================

def run_tsignal(input_fasta, output_csv, verbose=False):
    """Run TSignal prediction using Singularity container."""

    # Determine container runtime
    runtime = None
    for cmd in ["apptainer", "singularity"]:
        try:
            subprocess.run([cmd, "--version"], capture_output=True)
            runtime = cmd
            break
        except FileNotFoundError:
            continue

    if not runtime:
        raise RuntimeError("Neither 'apptainer' nor 'singularity' found")

    # Create temporary working directory
    work_dir = Path(input_fasta).parent.resolve()
    input_name = Path(input_fasta).name
    output_name = Path(output_csv).name

    # Copy input to sp_data location expected by TSignal
    sp_data_dir = SOURCE_DIR / "sp_data"
    sp_data_dir.mkdir(exist_ok=True)

    # Copy model if not already there
    model_dest = sp_data_dir / MODEL_NAME
    if not model_dest.exists():
        log(f"Linking model to sp_data/")
        os.symlink(MODEL_PATH, model_dest)

    # Copy input fasta
    import shutil
    input_dest = sp_data_dir / input_name
    shutil.copy(input_fasta, input_dest)

    # Build command
    cmd = [
        runtime, "run",
        "--bind", f"{SOURCE_DIR}:/app/TSignal",
        str(CONTAINER_PATH),
        "--test_seqs", f"sp_data/{input_name}",
        "--test_mdl", MODEL_NAME,
        "--tune_bert",
        "--train_only_decoder",
        "--output_file", f"sp_data/{output_name}"
    ]

    if verbose:
        cmd.append("--verbouse")

    log(f"Running TSignal on {input_name}...")
    log(f"Command: {' '.join(cmd)}")

    # Run TSignal
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(SOURCE_DIR)
    )

    if result.returncode != 0:
        log(f"TSignal stderr: {result.stderr}", "ERROR")
        raise RuntimeError(f"TSignal failed with return code {result.returncode}")

    # Move output to desired location
    output_src = sp_data_dir / output_name
    if output_src.exists():
        shutil.move(str(output_src), output_csv)
        log(f"Predictions saved to: {output_csv}")
    else:
        raise RuntimeError(f"TSignal did not produce output file: {output_src}")

    # Cleanup temporary input
    input_dest.unlink()

    return output_csv


# ============================================================================
# Signal Peptide Processing
# ============================================================================

def parse_predictions(predictions_csv):
    """Parse TSignal predictions CSV file."""
    predictions = {}

    with open(predictions_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row['seqs'][:60]  # TSignal uses first 60aa
            predictions[seq] = {
                'pred_lbls': row['pred_lbls'],
                'sp_type_prob': float(row['SP_type_prob']),
                'cs_prob': row.get('CS_prob', 'N/A')
            }

    return predictions


def find_cleavage_site(pred_labels):
    """
    Find cleavage site position from prediction labels.
    Returns 1-based position of the last SP residue.
    """
    sp_char = pred_labels[0]

    # Only process if it's a signal peptide type
    if sp_char not in ['S', 'T', 'L', 'W', 'P']:
        return 0

    cs_pos = 0
    for i, c in enumerate(pred_labels):
        if c == sp_char:
            cs_pos = i + 1  # 1-based position
        elif c in ['O', 'I', 'M', 'E']:
            break

    return cs_pos


def process_sequences(sequences, predictions, threshold=0.9):
    """
    Process sequences using predictions.

    Returns:
        annotations: List of annotation dicts
        processed: Dict of {header: processed_sequence}
        signal_peptides: Dict of {header: sp_sequence}
    """
    annotations = []
    processed = {}
    signal_peptides = {}

    for header, full_seq in sequences.items():
        first60 = full_seq[:60]

        if first60 in predictions:
            pred = predictions[first60]
            pred_labels = pred['pred_lbls']
            sp_prob = pred['sp_type_prob']
            cs_prob = pred['cs_prob']

            # Determine SP type
            sp_char = pred_labels[0]
            sp_type = SP_TYPES.get(sp_char, 'Unknown')

            # Find cleavage site
            cs_pos = find_cleavage_site(pred_labels)

            # Decide whether to remove SP
            is_sp = sp_char in ['S', 'T', 'L', 'W', 'P']
            high_confidence = sp_prob > threshold
            has_cleavage = cs_pos > 0

            if is_sp and high_confidence and has_cleavage:
                processed_seq = full_seq[cs_pos:]
                sp_sequence = full_seq[:cs_pos]
                sp_removed = True
            else:
                processed_seq = full_seq
                sp_sequence = ''
                sp_removed = False

            annotations.append({
                'sequence_id': header,
                'sp_type': sp_type,
                'sp_probability': sp_prob,
                'cleavage_site': cs_pos if sp_removed else 'N/A',
                'cs_probability': cs_prob,
                'sp_removed': sp_removed,
                'sp_sequence': sp_sequence,
                'original_length': len(full_seq),
                'mature_length': len(processed_seq),
                'prediction_labels': pred_labels[:30] + '...' if len(pred_labels) > 30 else pred_labels
            })

            processed[header] = processed_seq
            if sp_removed:
                signal_peptides[header] = sp_sequence

        else:
            # No prediction found (different N-terminus)
            annotations.append({
                'sequence_id': header,
                'sp_type': 'NO_PREDICTION',
                'sp_probability': 'N/A',
                'cleavage_site': 'N/A',
                'cs_probability': 'N/A',
                'sp_removed': False,
                'sp_sequence': '',
                'original_length': len(full_seq),
                'mature_length': len(full_seq),
                'prediction_labels': 'N/A'
            })
            processed[header] = full_seq

    return annotations, processed, signal_peptides


def write_annotations(annotations, output_path):
    """Write annotations to TSV file."""
    fieldnames = [
        'sequence_id', 'sp_type', 'sp_probability', 'cleavage_site',
        'cs_probability', 'sp_removed', 'sp_sequence', 'original_length',
        'mature_length', 'prediction_labels'
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(annotations)


def write_summary(annotations, output_path, input_file, threshold):
    """Write processing summary."""
    total = len(annotations)
    sp_removed = sum(1 for a in annotations if a['sp_removed'])
    no_prediction = sum(1 for a in annotations if a['sp_type'] == 'NO_PREDICTION')

    # Count by SP type
    sp_counts = defaultdict(int)
    for a in annotations:
        sp_counts[a['sp_type']] += 1

    with open(output_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("TSignal Signal Peptide Processing Summary\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Input file: {input_file}\n")
        f.write(f"Probability threshold: {threshold}\n")
        f.write(f"Processing date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("-" * 40 + "\n")
        f.write("Results\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total sequences: {total}\n")
        f.write(f"Signal peptides removed: {sp_removed}\n")
        f.write(f"Sequences unchanged: {total - sp_removed}\n")
        f.write(f"  - Below threshold: {total - sp_removed - no_prediction}\n")
        f.write(f"  - No prediction: {no_prediction}\n\n")

        f.write("-" * 40 + "\n")
        f.write("Signal Peptide Types Detected\n")
        f.write("-" * 40 + "\n")
        for sp_type, count in sorted(sp_counts.items(), key=lambda x: -x[1]):
            pct = count / total * 100
            f.write(f"  {sp_type}: {count} ({pct:.1f}%)\n")

        f.write("\n" + "=" * 60 + "\n")

    return sp_removed, total


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="TSignal batch processing: predict and remove signal peptides",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i sequences.fasta -o output/
  %(prog)s -i sequences.fasta -o output/ --threshold 0.8
  %(prog)s --check

Output files:
  tsignal_predictions.csv      - Raw TSignal predictions
  sequences_sp_removed.fasta   - Sequences with SPs removed
  signal_peptides_only.fasta   - Extracted SP sequences
  signal_peptide_annotations.tsv - Detailed annotations
  processing_summary.txt       - Summary statistics
        """
    )

    parser.add_argument('-i', '--input', type=str,
                        help='Input FASTA file')
    parser.add_argument('-o', '--output-dir', type=str, default='output',
                        help='Output directory (default: output)')
    parser.add_argument('--threshold', type=float, default=0.9,
                        help='SP probability threshold for removal (default: 0.9)')
    parser.add_argument('--skip-prediction', action='store_true',
                        help='Skip TSignal prediction (use existing predictions.csv)')
    parser.add_argument('--predictions', type=str,
                        help='Path to existing predictions CSV (with --skip-prediction)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Verbose output')
    parser.add_argument('--check', action='store_true',
                        help='Check installation and exit')

    args = parser.parse_args()

    # Check installation
    if args.check:
        success = check_installation()
        sys.exit(0 if success else 1)

    # Validate arguments
    if not args.input:
        parser.error("--input is required")

    if not os.path.exists(args.input):
        log(f"Input file not found: {args.input}", "ERROR")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define output paths
    predictions_csv = output_dir / "tsignal_predictions.csv"
    sequences_out = output_dir / "sequences_sp_removed.fasta"
    sp_only_out = output_dir / "signal_peptides_only.fasta"
    annotations_out = output_dir / "signal_peptide_annotations.tsv"
    summary_out = output_dir / "processing_summary.txt"

    # Step 1: Read input sequences
    log(f"Reading input: {args.input}")
    sequences = read_fasta(args.input)
    log(f"Loaded {len(sequences)} sequences")

    # Step 2: Run TSignal prediction
    if args.skip_prediction and args.predictions:
        log(f"Using existing predictions: {args.predictions}")
        predictions_csv = Path(args.predictions)
    else:
        if not check_installation():
            sys.exit(1)
        run_tsignal(args.input, str(predictions_csv), verbose=args.verbose)

    # Step 3: Parse predictions
    log("Parsing predictions...")
    predictions = parse_predictions(predictions_csv)
    log(f"Loaded {len(predictions)} unique predictions")

    # Step 4: Process sequences
    log(f"Processing sequences (threshold={args.threshold})...")
    annotations, processed, signal_peptides = process_sequences(
        sequences, predictions, threshold=args.threshold
    )

    # Step 5: Write outputs
    log("Writing output files...")

    write_fasta(processed, sequences_out)
    log(f"  -> {sequences_out}")

    if signal_peptides:
        write_fasta(signal_peptides, sp_only_out)
        log(f"  -> {sp_only_out}")

    write_annotations(annotations, annotations_out)
    log(f"  -> {annotations_out}")

    sp_removed, total = write_summary(
        annotations, summary_out, args.input, args.threshold
    )
    log(f"  -> {summary_out}")

    # Final summary
    log("=" * 50)
    log(f"COMPLETE: {sp_removed}/{total} signal peptides removed")
    log("=" * 50)


if __name__ == "__main__":
    main()
