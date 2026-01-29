# TSignal Toolkit

**A ready-to-use toolkit for predicting and removing signal peptides from protein sequences.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## About This Toolkit

This toolkit is a **user-friendly repackaging** of [TSignal](https://github.com/Dumitrescu-Alexandru/TSignal), a state-of-the-art transformer-based signal peptide predictor developed by Dumitrescu et al. (2022). The original TSignal software requires manual setup of dependencies, conda environments, and careful configuration. This toolkit wraps everything into a portable Singularity container with simple batch processing scripts, making it easy to:

1. **Predict signal peptides** in your protein sequences
2. **Remove signal peptides** to obtain mature protein sequences
3. **Extract signal peptides** for separate analysis
4. **Generate detailed annotations** for every sequence

**If you use this toolkit, please cite the original TSignal publication** (see [Citation](#citation)).

---

## Background: What Are Signal Peptides?

**Signal peptides (SPs)** are short amino acid sequences (typically 15-30 residues) located at the N-terminus of proteins destined for secretion or membrane insertion. They act as "address labels" that direct newly synthesized proteins to the secretory pathway.

### The Secretion Process

1. A ribosome begins translating an mRNA encoding a secretory protein
2. The emerging signal peptide is recognized by the Signal Recognition Particle (SRP)
3. The ribosome-SRP complex docks at the ER membrane (eukaryotes) or plasma membrane (prokaryotes)
4. The protein is threaded through a translocon channel into the ER lumen or periplasm
5. **Signal peptidases cleave off the signal peptide**, releasing the mature protein

### Why Remove Signal Peptides for Analysis?

Since signal peptides are **cleaved and degraded** during protein maturation, the functional protein that exists in the cell lacks this N-terminal sequence. For many bioinformatics analyses, working with the **mature protein sequence** (without SP) is preferable:

- **Phylogenetic analysis**: Signal peptides evolve differently than the mature protein and can introduce noise into alignments and tree reconstruction
- **Structure prediction**: AlphaFold and similar tools work better on mature sequences
- **Domain analysis**: Signal peptides can interfere with domain detection
- **Sequence alignments**: Mature regions align more reliably across homologs

### Signal Peptide Types

There are five major types of signal peptides, distinguished by their pathway and the peptidase that cleaves them:

| Type | Pathway | Peptidase | Description |
|------|---------|-----------|-------------|
| **Sec/SPase I** | Sec | Signal Peptidase I | Most common type; secreted proteins and many membrane proteins |
| **Sec/SPase II** | Sec | Signal Peptidase II | Lipoproteins; have a conserved cysteine that is lipid-modified |
| **Tat/SPase I** | Tat | Signal Peptidase I | Twin-arginine translocation; proteins that fold before transport |
| **Tat/SPase II** | Tat | Signal Peptidase II | Tat-transported lipoproteins |
| **Sec/SPase IV** | Sec | Prepilin Peptidase | Type IV pilins and related proteins (bacteria) |

TSignal can predict all five types with high accuracy.

---

## About TSignal

TSignal is a deep learning model for signal peptide prediction, developed by Dumitrescu et al. and described in their 2022 preprint:

> **Dumitrescu, A., Sharan, M., Gîlcă, G., Henkel, M., Bernhofer, M., Rost, B., & Achim, A. (2022).** TSignal: A transformer model for signal peptide prediction. *bioRxiv*. https://doi.org/10.1101/2022.06.02.493958

### Key Features of TSignal

- **Transformer architecture**: Adapts the encoder-decoder architecture from "Attention is All You Need" (Vaswani et al., 2017) for sequence labeling
- **Transfer learning**: Uses pre-trained ProtBERT weights for protein language understanding
- **All five SP types**: Predicts Sec/SPase I, Sec/SPase II, Tat/SPase I, Tat/SPase II, and Sec/SPase IV
- **Cleavage site prediction**: Identifies the exact position where the signal peptide is cleaved
- **High accuracy**: Competitive with or outperforms SignalP 6.0 on benchmark datasets

### Original TSignal Repository

- **GitHub**: https://github.com/Dumitrescu-Alexandru/TSignal
- **Preprint**: https://www.biorxiv.org/content/10.1101/2022.06.02.493958

---

## What This Toolkit Adds

This toolkit repackages TSignal with:

1. **Singularity container**: All dependencies pre-configured for portability across systems
2. **Batch processing scripts**: Simple command-line interface for processing FASTA files
3. **Automatic SP removal**: Parses predictions and generates mature protein sequences
4. **Detailed output**: TSV annotations, extracted SPs, and summary statistics
5. **HPC support**: SLURM job submission scripts included

**No changes were made to the core TSignal prediction algorithm.** This toolkit simply makes it easier to use.

---

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Detailed Usage](#detailed-usage)
4. [Output Files](#output-files)
5. [Advanced Options](#advanced-options)
6. [Running on HPC Clusters](#running-on-hpc-clusters)
7. [Troubleshooting](#troubleshooting)
8. [Citation](#citation)
9. [License](#license)

---

## Installation

### Step 1: Install Singularity/Apptainer

**Singularity or Apptainer** (container runtime) is required.

Check if you have it:
```bash
singularity --version   # or
apptainer --version
```

If not installed:
- **Linux (Ubuntu/Debian):** `sudo apt install apptainer`
- **Linux (CentOS/RHEL):** `sudo yum install apptainer`
- **HPC clusters:** Usually pre-installed - try `module load singularity` or `module load apptainer`
- **Mac:** Use a Linux VM or Docker with Singularity

**Python 3.6+** is also required:
```bash
python3 --version
```

### Step 2: Download the Toolkit

**Option A: Download Release (Recommended)**

Download the latest release which includes the Singularity container:

```bash
# Download and extract the release (includes container)
wget https://github.com/invertome/TSignal_toolkit/releases/download/v1.0.0/TSignal_toolkit-v1.0.0.tar.gz
tar -xzf TSignal_toolkit-v1.0.0.tar.gz
cd TSignal_toolkit-v1.0.0

# Download only the pre-trained model (~1.9 GB)
./scripts/download_models.sh
```

**Option B: Clone from Git**

```bash
git clone https://github.com/invertome/TSignal_toolkit.git
cd TSignal_toolkit

# Download container and model (~2.5 GB total)
./scripts/download_models.sh
```

### Step 3: Verify Installation

```bash
./scripts/tsignal_batch.py --check
```

Expected output:
```
[INFO] Container found: container/tsignal.sif (557.0 MB)
[INFO] Model found: model/deployment_*.pth (1.95 GB)
[INFO] All required files present!
[INFO] Installation check PASSED
```

---

## Quick Start

```bash
# Run on your sequences (does prediction + SP removal)
./scripts/tsignal_batch.py -i your_sequences.fasta -o results/

# Check results
ls results/
#   sequences_sp_removed.fasta      <- Mature proteins (SPs removed)
#   signal_peptides_only.fasta      <- Extracted SPs
#   signal_peptide_annotations.tsv  <- Detailed per-sequence info
#   processing_summary.txt          <- Summary statistics
```

---

## Detailed Usage

### Basic Command

```bash
./scripts/tsignal_batch.py -i INPUT.fasta -o OUTPUT_DIR/
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `-i, --input` | Input FASTA file (required) | - |
| `-o, --output-dir` | Output directory | `output/` |
| `--threshold` | Probability threshold for SP removal | `0.9` |
| `-t, --threads` | Number of CPU threads | all available |
| `--batch-size` | Process in batches (for large files) | disabled |
| `-v, --verbose` | Show detailed progress | off |
| `--check` | Verify installation | - |
| `--version` | Show version number | - |

### Examples

```bash
# Basic usage
./scripts/tsignal_batch.py -i proteins.fasta -o results/

# Lower threshold (more permissive, removes more SPs)
./scripts/tsignal_batch.py -i proteins.fasta -o results/ --threshold 0.8

# Higher threshold (more conservative, only very confident predictions)
./scripts/tsignal_batch.py -i proteins.fasta -o results/ --threshold 0.95

# Limit to 4 CPU threads
./scripts/tsignal_batch.py -i proteins.fasta -o results/ --threads 4

# Process large file in batches of 5000 sequences
./scripts/tsignal_batch.py -i large_dataset.fasta -o results/ --batch-size 5000

# Verbose mode
./scripts/tsignal_batch.py -i proteins.fasta -o results/ --verbose
```

### Step-by-Step Processing

For more control, run prediction and processing separately:

```bash
# Step 1: Run TSignal prediction only
./scripts/tsignal_predict.sh sequences.fasta predictions.csv

# Step 2: Process predictions (can try different thresholds)
./scripts/process_predictions.py \
    --predictions predictions.csv \
    --fasta sequences.fasta \
    --output-dir results/ \
    --threshold 0.9
```

---

## Output Files

### `sequences_sp_removed.fasta`

Mature protein sequences with signal peptides removed. Sequences below the confidence threshold are included unchanged.

### `signal_peptides_only.fasta`

Only the extracted signal peptide sequences (for sequences where SP was removed).

### `signal_peptide_annotations.tsv`

Detailed tab-separated annotations:

| Column | Description |
|--------|-------------|
| `sequence_id` | FASTA header |
| `sp_type` | Signal peptide type |
| `sp_probability` | Confidence score (0-1) |
| `cleavage_site` | Cleavage position |
| `sp_removed` | Whether SP was removed |
| `sp_sequence` | The signal peptide sequence |
| `original_length` | Original sequence length |
| `mature_length` | Length after SP removal |

### `processing_summary.txt`

Human-readable summary with statistics.

---

## Advanced Options

### Processing Existing Predictions

Re-process with different thresholds without re-running TSignal:

```bash
./scripts/process_predictions.py \
    -p results/tsignal_predictions.csv \
    -f original.fasta \
    -o results_0.8/ \
    --threshold 0.8
```

### Batch Processing

```bash
for f in data/*.fasta; do
    name=$(basename "$f" .fasta)
    ./scripts/tsignal_batch.py -i "$f" -o "results/${name}/"
done
```

---

## Running on HPC Clusters

### SLURM Submission

```bash
# Edit settings
nano scripts/tsignal_slurm.sbatch

# Submit
sbatch scripts/tsignal_slurm.sbatch
```

### Key SLURM Settings

```bash
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00

INPUT_FASTA="your_sequences.fasta"
OUTPUT_DIR="output"
THRESHOLD="0.9"
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Container not found | Run `./scripts/download_models.sh` |
| Model not found | Download from Dropbox (see Installation) |
| singularity not found | `module load singularity` or install Apptainer |
| Permission denied | `chmod +x scripts/*.py scripts/*.sh` |
| Slow performance | Use SLURM script on HPC cluster |

---

## Citation

**If you use this toolkit, please cite the original TSignal publication:**

```bibtex
@article{dumitrescu2022tsignal,
  title={TSignal: A transformer model for signal peptide prediction},
  author={Dumitrescu, Alexandru and Sharan, Moksh and G\^{i}lc\u{a}, Andrei and
          Henkel, Maria and Bernhofer, Michael and Rost, Burkhard and Achim, Alin},
  journal={bioRxiv},
  year={2022},
  doi={10.1101/2022.06.02.493958}
}
```

**Also cite SignalP 6.0 for the training data:**

```bibtex
@article{teufel2022signalp,
  title={SignalP 6.0 predicts all five types of signal peptides using protein language models},
  author={Teufel, Felix and Almagro Armenteros, Jos{\'e} Juan and Johansen, Alexander Rosenberg and
          G{\'i}slason, Magn{\'u}s Halld{\'o}r and Piber, Silas Irber and Tsirigos, Konstantinos D and
          Winther, Ole and Brunak, S{\o}ren and von Heijne, Gunnar and Nielsen, Henrik},
  journal={Nature Biotechnology},
  volume={40},
  number={7},
  pages={1023--1025},
  year={2022},
  publisher={Nature Publishing Group}
}
```

---

## License

This toolkit is released under the **MIT License**.

The original TSignal software by Dumitrescu et al. is also released under the MIT License.

See [LICENSE](LICENSE) for details.

---

## References

- **TSignal GitHub**: https://github.com/Dumitrescu-Alexandru/TSignal
- **TSignal Preprint**: https://www.biorxiv.org/content/10.1101/2022.06.02.493958
- **SignalP 6.0**: https://services.healthtech.dtu.dk/service.php?SignalP-6.0
- **ProtTrans**: https://github.com/agemagician/ProtTrans

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

---

*Toolkit assembled by Jorge L. Perez-Moreno, Ph.D. (jperezmoreno@umass.edu)*
