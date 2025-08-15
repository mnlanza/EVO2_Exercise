# EVO2 Analysis Pipeline

This repository contains scripts for analyzing gene variants using the EVO2 model. The pipeline consists of three main analyses:
1. Standard variant analysis
2. Shifting positions analysis
3. Reverse complement analysis

## Directory Structure
```
.
├── scripts/
│   ├── Ex1_main.R               # Main R script for analysis
│   ├── evo2_analysis_functions.R # Helper functions for plots
│   ├── handle_EVO2_output.r     # Functions to process EVO2 output
│   ├── shift_sequences.py       # Script to create shifted sequences
│   ├── generate_r_cods_v1.py    # Generate reverse complement sequences
│   └── frameshift_fasta.py      # Script to handle frameshifts
├── raw_data/
│   ├── gene_variants.fasta      # Original gene variants
│   ├── codon83_variants.fasta   # Input variants for reverse complement
│   ├── shifted_variants.fasta   # Shifted referece sequence variants shifted 1 & 2 nucleotides
│   ├── reverse_complements_83.fasta # Generated reverse complement variants for all codons at position 83
│   └── input_*_logits.npy      # EVO2 logits files
├── plots/
│   ├── ex1_plots/              # Standard analysis plots
│   ├── ex1_shifting_plots/     # Shifting analysis plots
│   └── rc_Ex1/                 # Reverse complement analysis plots
└── run_analysis.sh             # Main runner script
```

## Workflow

### 1. Prepare FASTA Files

#### Standard Variants
- Using provided `gene_variants.fasta` containing variants:
  - 83_S1 (reference)
  - 83_L
  - 83_S2
  - 83_Stop

#### Shifted Variants
`raw_data/shifted_variants.fasta` is provided. If you need to regenerate shifted sequences, you can adapt `scripts/frameshift_fasta.py`. This analysis didn't prove much though as the results were as expected: the extra nucleotides before the start codon 1,2 didn't effect the result. 

#### Reverse Complement Variants
Generate reverse complement sequences using:
```bash
cd raw_data
python ../scripts/generate_r_cods_v1.py
```
This takes `codon83_variants.fasta` as input and creates `reverse_complements_83.fasta` with reverse complement sequences:
- Each sequence ID gets "_rc" appended to indicate it is a reverse complement.

### 2. Run EVO2 Model
For each FASTA file, run the EVO2 model to generate logits files in the `raw_data` directory:
- Standard variants: input_83_S1_logits.npy, input_83_L_logits.npy, etc.
- Shifted variants: input_83_S1_shifted1r_logits.npy, input_83_S1_shifted2r_logits.npy, etc.
- Reverse complement: input_83_AGC_rc_logits.npy, input_83_ATC_rc_logits.npy, etc.

### 3. Run Analysis
Execute all analyses using:
```bash
./run_analysis.sh
```

This script runs three analyses:
1. Standard variant analysis
2. Shifting positions analysis (with reverse complement flag)
3. Reverse complement analysis

You can also run the main R script directly:
```bash
Rscript scripts/Ex1_main.R \
  --fasta raw_data/gene_variants.fasta \
  --variants 83_S1,83_L,83_S2,83_Stop \
  --input_logits_dir raw_data \
  --output_dir plots/ex1_plots \
  --highlight 247,248,249 \
  --index_start 200 --index_end 300 \
  --reverse_comp \
  --only_single
```

CLI options:
- `--fasta`: path to FASTA file
- `--variants`: comma-separated IDs in FASTA order (first is the reference)
- `--input_logits_dir`: directory of `input_*_logits.npy`
- `--output_dir`: where plots are written
- `--highlight`: comma-separated positions for the mutated region
- `--index_start` / `--index_end`: plotting range
- `--reverse_comp`: flip coordinates for reverse complements
- `--only_single`: generate only single-variant plots (difference plots are still saved in stacked form)

## Output
Results are organized in three directories:
- `plots/ex1_plots/`: Standard variant analysis
- `plots/ex1_shifting_plots/`: Shifting positions analysis
- `plots/rc_Ex1/`: Reverse complement analysis

Each directory contains:
- `single/`: Stacked single-variant plots (Entropy over −Log-likelihood) per variant, one PDF per variant
- `compare/`: Stacked difference plots across variants (two PDFs):
  - `entropy_diffs_stacked.pdf`
  - `ll_diffs_stacked.pdf`
- `tot_ll/`: Log-likelihood heatmaps

Notes on plots:
- Differences are computed as `target − reference`, where the reference is the first variant in `--variants`.
- Log-likelihood is displayed as `−Log-likelihood` (capital L) throughout.
- The highlighted vertical rectangles are labeled as `Mutated region` in the legend.
- In stacked figures, only the bottom panel shows the x‑axis to reduce clutter.

## Dependencies
- R packages:
  - optparse
  - Biostrings
  - ggplot2
  - patchwork
  - reticulate (for Python integration)
- Python packages:
  - Biopython
  - numpy
  - argparse
- EVO2 model and dependencies
