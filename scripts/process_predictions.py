#!/usr/bin/env python3
"""
Process TSignal Predictions
============================
Parse TSignal output and remove signal peptides from sequences.

This script can be used standalone when you already have TSignal predictions.

Usage:
    ./process_predictions.py --predictions predictions.csv --fasta sequences.fasta --output-dir output/

Author: Jorge L. Perez-Moreno, Ph.D. (jperezmoreno@umass.edu)
Based on TSignal by Dumitrescu et al. (2022)
https://github.com/Dumitrescu-Alexandru/TSignal
"""

__version__ = "1.0.0"

import argparse
import csv
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path


# Signal peptide type mapping
SP_TYPES = {
    'S': 'Sec/SPase_I',    # Standard secretory pathway
    'T': 'Tat/SPase_I',    # Twin-arginine translocation
    'L': 'Sec/SPase_II',   # Lipoprotein signal peptide
    'W': 'Tat/SPase_II',   # Tat lipoprotein
    'P': 'Sec/SPase_IV',   # Prepilin-like
    'I': 'No_SP',          # No signal peptide
    'O': 'Other',          # Post-cleavage region
    'M': 'Membrane',       # Transmembrane
    'E': 'End'             # Sequence end
}


def log(msg, level="INFO"):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {level}: {msg}", file=sys.stderr)


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
                current_header = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line)

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


def parse_predictions(predictions_csv):
    """
    Parse TSignal predictions CSV file.

    Returns dict mapping first 60aa to prediction info.
    """
    predictions = {}

    with open(predictions_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # TSignal uses first 60aa for prediction
            seq = row['seqs'][:60]
            predictions[seq] = {
                'full_seq': row['seqs'],
                'pred_lbls': row['pred_lbls'],
                'sp_type_prob': float(row['SP_type_prob']),
                'cs_prob': row.get('CS_prob', 'N/A')
            }

    return predictions


def find_cleavage_site(pred_labels):
    """
    Find cleavage site position from prediction labels.

    TSignal labels each residue with:
    - S/T/L/W/P = Signal peptide residue (different types)
    - O = Other (mature protein)
    - I = Cytoplasmic (no SP)
    - M = Membrane

    Returns 1-based position of the last SP residue (cleavage occurs after this).
    """
    sp_char = pred_labels[0]

    # Only look for cleavage site if it's a signal peptide type
    if sp_char not in ['S', 'T', 'L', 'W', 'P']:
        return 0

    cs_pos = 0
    for i, c in enumerate(pred_labels):
        if c == sp_char:
            cs_pos = i + 1  # 1-based position
        elif c in ['O', 'I', 'M', 'E']:
            # End of signal peptide region
            break

    return cs_pos


def process_sequences(sequences, predictions, threshold=0.9):
    """
    Process sequences using TSignal predictions.

    Args:
        sequences: Dict of {header: full_sequence}
        predictions: Dict of {first60aa: prediction_dict}
        threshold: Probability threshold for SP removal (default 0.9)

    Returns:
        annotations: List of annotation dicts for each sequence
        processed: Dict of {header: processed_sequence} (SPs removed)
        signal_peptides: Dict of {header: sp_sequence} (extracted SPs)
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

            # Determine SP type from first prediction character
            sp_char = pred_labels[0]
            sp_type = SP_TYPES.get(sp_char, 'Unknown')

            # Find cleavage site
            cs_pos = find_cleavage_site(pred_labels)

            # Criteria for SP removal:
            # 1. Is a signal peptide type (S, T, L, W, or P)
            # 2. High confidence (probability > threshold)
            # 3. Has valid cleavage site
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
                'sp_probability': round(sp_prob, 4),
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
            # No prediction found (sequence not in predictions file)
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
    """Write detailed annotations to TSV file."""
    fieldnames = [
        'sequence_id', 'sp_type', 'sp_probability', 'cleavage_site',
        'cs_probability', 'sp_removed', 'sp_sequence', 'original_length',
        'mature_length', 'prediction_labels'
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(annotations)


def print_summary(annotations, threshold):
    """Print processing summary to stderr."""
    total = len(annotations)
    sp_removed = sum(1 for a in annotations if a['sp_removed'])
    no_prediction = sum(1 for a in annotations if a['sp_type'] == 'NO_PREDICTION')

    # Count by SP type
    sp_counts = defaultdict(int)
    for a in annotations:
        sp_counts[a['sp_type']] += 1

    log("=" * 50)
    log("PROCESSING SUMMARY")
    log("=" * 50)
    log(f"Probability threshold: {threshold}")
    log(f"Total sequences: {total}")
    log(f"Signal peptides removed: {sp_removed} ({sp_removed/total*100:.1f}%)")
    log(f"Sequences unchanged: {total - sp_removed}")
    log("")
    log("By SP type:")
    for sp_type, count in sorted(sp_counts.items(), key=lambda x: -x[1]):
        log(f"  {sp_type}: {count}")
    log("=" * 50)


def main():
    parser = argparse.ArgumentParser(
        description="Process TSignal predictions and remove signal peptides",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Signal Peptide Types:
  Sec/SPase_I  - Standard secretory pathway (most common)
  Sec/SPase_II - Lipoprotein signal peptide
  Tat/SPase_I  - Twin-arginine translocation pathway
  Tat/SPase_II - Tat lipoprotein
  Sec/SPase_IV - Prepilin-like (bacterial)
  No_SP        - No signal peptide detected

Output Files:
  sequences_sp_removed.fasta    - Sequences with SPs removed
  signal_peptides_only.fasta    - Extracted SP sequences only
  signal_peptide_annotations.tsv - Detailed annotations
        """
    )

    parser.add_argument('--predictions', '-p', required=True,
                        help='TSignal predictions CSV file')
    parser.add_argument('--fasta', '-f', required=True,
                        help='Original input FASTA file')
    parser.add_argument('--output-dir', '-o', default='output',
                        help='Output directory (default: output)')
    parser.add_argument('--threshold', '-t', type=float, default=0.9,
                        help='SP probability threshold for removal (default: 0.9)')
    parser.add_argument('--prefix', default='',
                        help='Prefix for output filenames')

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.predictions).exists():
        log(f"Predictions file not found: {args.predictions}", "ERROR")
        sys.exit(1)

    if not Path(args.fasta).exists():
        log(f"FASTA file not found: {args.fasta}", "ERROR")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define output paths
    prefix = f"{args.prefix}_" if args.prefix else ""
    sequences_out = output_dir / f"{prefix}sequences_sp_removed.fasta"
    sp_only_out = output_dir / f"{prefix}signal_peptides_only.fasta"
    annotations_out = output_dir / f"{prefix}signal_peptide_annotations.tsv"

    # Step 1: Read input sequences
    log(f"Reading sequences: {args.fasta}")
    sequences = read_fasta(args.fasta)
    log(f"Loaded {len(sequences)} sequences")

    # Step 2: Parse predictions
    log(f"Parsing predictions: {args.predictions}")
    predictions = parse_predictions(args.predictions)
    log(f"Loaded {len(predictions)} predictions")

    # Step 3: Process sequences
    log(f"Processing (threshold={args.threshold})...")
    annotations, processed, signal_peptides = process_sequences(
        sequences, predictions, threshold=args.threshold
    )

    # Step 4: Write outputs
    log("Writing output files...")

    write_fasta(processed, sequences_out)
    log(f"  -> {sequences_out}")

    if signal_peptides:
        write_fasta(signal_peptides, sp_only_out)
        log(f"  -> {sp_only_out} ({len(signal_peptides)} sequences)")
    else:
        log("  -> No signal peptides to extract")

    write_annotations(annotations, annotations_out)
    log(f"  -> {annotations_out}")

    # Print summary
    print_summary(annotations, args.threshold)


if __name__ == "__main__":
    main()
